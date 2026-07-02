#include <GHTTPConnection.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <random>
#include <sstream>

#include <netdb.h>
#include <sys/socket.h>
#include <unistd.h>

GHTTPConnection::GHTTPConnection()
  : fSocket(-1), fStopRequested(false), fConnected(false) {
}

GHTTPConnection::~GHTTPConnection() {
  Stop();
}

bool GHTTPConnection::Start(const std::string& url) {
  fUrl = url;
  if(!ParseUrl(url, fParsedUrl)) {
    std::cout << "Could not parse live histogram URL: " << url << std::endl;
    return false;
  }

  if(!ConnectSocket()) {
    std::cout << "Could not connect to live histogram server "
              << fParsedUrl.host << ":" << fParsedUrl.port << std::endl;
    return false;
  }

  if(!SendHandshake()) {
    std::cout << "Live histogram websocket handshake failed for "
              << fUrl << std::endl;
    CloseSocket();
    return false;
  }

  fStopRequested = false;
  fConnected = true;

  std::cout << "\tconnected live histogram source: " << fUrl << std::endl
            << "\twebsocket endpoint: ws://" << fParsedUrl.host
            << ":" << fParsedUrl.port << fParsedUrl.path << std::endl;

  SendText("LIST");
  SendText("SUBSCRIBE *");

  fThread = std::thread(&GHTTPConnection::ReceiverLoop, this);
  return true;
}

void GHTTPConnection::Stop() {
  fStopRequested = true;
  CloseSocket();

  if(fThread.joinable()) {
    fThread.join();
  }

  fConnected = false;
}

bool GHTTPConnection::ParseUrl(const std::string& url, ParsedUrl& parsed) {
  std::string work = url;
  const std::string http_prefix = "http://";
  const std::string ws_prefix = "ws://";

  if(work.compare(0, http_prefix.size(), http_prefix) == 0) {
    work = work.substr(http_prefix.size());
  } else if(work.compare(0, ws_prefix.size(), ws_prefix) == 0) {
    work = work.substr(ws_prefix.size());
  } else {
    return false;
  }

  std::string::size_type slash = work.find('/');
  std::string host_port = slash == std::string::npos ? work : work.substr(0, slash);
  parsed.path = slash == std::string::npos ? "/live/" : work.substr(slash);

  if(parsed.path.empty() || parsed.path == "/") {
    parsed.path = "/live/";
  }
  if(parsed.path.back() != '/') {
    parsed.path += "/";
  }
  if(parsed.path.find("root.websocket") == std::string::npos) {
    parsed.path += "root.websocket";
  }

  std::string::size_type colon = host_port.rfind(':');
  if(colon == std::string::npos) {
    parsed.host = host_port;
    parsed.port = "80";
  } else {
    parsed.host = host_port.substr(0, colon);
    parsed.port = host_port.substr(colon + 1);
  }

  return !parsed.host.empty() && !parsed.port.empty();
}

std::string GHTTPConnection::Base64Encode(const unsigned char* data, size_t length) {
  static const char alphabet[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  std::string output;
  output.reserve(((length + 2) / 3) * 4);

  for(size_t i = 0; i < length; i += 3) {
    unsigned int value = static_cast<unsigned int>(data[i]) << 16;
    if(i + 1 < length) {
      value |= static_cast<unsigned int>(data[i + 1]) << 8;
    }
    if(i + 2 < length) {
      value |= static_cast<unsigned int>(data[i + 2]);
    }

    output.push_back(alphabet[(value >> 18) & 0x3f]);
    output.push_back(alphabet[(value >> 12) & 0x3f]);
    output.push_back(i + 1 < length ? alphabet[(value >> 6) & 0x3f] : '=');
    output.push_back(i + 2 < length ? alphabet[value & 0x3f] : '=');
  }

  return output;
}

std::string GHTTPConnection::MakeWebSocketKey() {
  std::array<unsigned char, 16> bytes;
  std::random_device random;
  for(unsigned char& byte : bytes) {
    byte = static_cast<unsigned char>(random());
  }
  return Base64Encode(bytes.data(), bytes.size());
}

bool GHTTPConnection::ConnectSocket() {
  struct addrinfo hints;
  std::memset(&hints, 0, sizeof(hints));
  hints.ai_family = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;

  struct addrinfo* result = nullptr;
  int status = getaddrinfo(fParsedUrl.host.c_str(),
                           fParsedUrl.port.c_str(),
                           &hints, &result);
  if(status != 0) {
    return false;
  }

  for(struct addrinfo* rp = result; rp; rp = rp->ai_next) {
    int candidate = socket(rp->ai_family, rp->ai_socktype, rp->ai_protocol);
    if(candidate < 0) {
      continue;
    }
    if(connect(candidate, rp->ai_addr, rp->ai_addrlen) == 0) {
      fSocket = candidate;
      break;
    }
    close(candidate);
  }

  freeaddrinfo(result);
  return fSocket >= 0;
}

bool GHTTPConnection::SendHandshake() {
  std::ostringstream request;
  request << "GET " << fParsedUrl.path << " HTTP/1.1\r\n"
          << "Host: " << fParsedUrl.host << ":" << fParsedUrl.port << "\r\n"
          << "Upgrade: websocket\r\n"
          << "Connection: Upgrade\r\n"
          << "Sec-WebSocket-Key: " << MakeWebSocketKey() << "\r\n"
          << "Sec-WebSocket-Version: 13\r\n\r\n";

  const std::string text = request.str();
  if(send(fSocket, text.data(), text.size(), 0) != static_cast<ssize_t>(text.size())) {
    return false;
  }

  std::string response;
  char buffer[1024];
  while(response.find("\r\n\r\n") == std::string::npos) {
    ssize_t nread = recv(fSocket, buffer, sizeof(buffer), 0);
    if(nread <= 0) {
      return false;
    }
    response.append(buffer, static_cast<size_t>(nread));
  }

  return response.find(" 101 ") != std::string::npos ||
         response.find(" 101\r\n") != std::string::npos;
}

bool GHTTPConnection::SendText(const std::string& text) {
  if(fSocket < 0) {
    return false;
  }

  std::vector<unsigned char> frame;
  frame.push_back(0x81);

  const size_t length = text.size();
  if(length < 126) {
    frame.push_back(static_cast<unsigned char>(0x80 | length));
  } else if(length <= 0xffff) {
    frame.push_back(0x80 | 126);
    frame.push_back(static_cast<unsigned char>((length >> 8) & 0xff));
    frame.push_back(static_cast<unsigned char>(length & 0xff));
  } else {
    frame.push_back(0x80 | 127);
    for(int shift = 56; shift >= 0; shift -= 8) {
      frame.push_back(static_cast<unsigned char>((length >> shift) & 0xff));
    }
  }

  const unsigned char mask[4] = {0x47, 0x52, 0x4f, 0x54};
  frame.insert(frame.end(), mask, mask + 4);
  for(size_t i = 0; i < text.size(); ++i) {
    frame.push_back(static_cast<unsigned char>(text[i]) ^ mask[i % 4]);
  }

  return send(fSocket, frame.data(), frame.size(), 0) == static_cast<ssize_t>(frame.size());
}

bool GHTTPConnection::ReadExact(void* data, size_t length) {
  char* out = static_cast<char*>(data);
  size_t done = 0;
  while(done < length && !fStopRequested) {
    ssize_t nread = recv(fSocket, out + done, length - done, 0);
    if(nread <= 0) {
      return false;
    }
    done += static_cast<size_t>(nread);
  }
  return done == length;
}

bool GHTTPConnection::ReadFrame(unsigned char& opcode, std::vector<char>& payload) {
  unsigned char header[2];
  if(!ReadExact(header, 2)) {
    return false;
  }

  opcode = header[0] & 0x0f;
  bool masked = (header[1] & 0x80) != 0;
  unsigned long long length = header[1] & 0x7f;
  if(length == 126) {
    unsigned char extended[2];
    if(!ReadExact(extended, 2)) {
      return false;
    }
    length = (static_cast<unsigned long long>(extended[0]) << 8) | extended[1];
  } else if(length == 127) {
    unsigned char extended[8];
    if(!ReadExact(extended, 8)) {
      return false;
    }
    length = 0;
    for(unsigned char c : extended) {
      length = (length << 8) | c;
    }
  }

  unsigned char mask[4] = {0, 0, 0, 0};
  if(masked && !ReadExact(mask, 4)) {
    return false;
  }

  payload.resize(static_cast<size_t>(length));
  if(length && !ReadExact(payload.data(), static_cast<size_t>(length))) {
    return false;
  }

  if(masked) {
    for(size_t i = 0; i < payload.size(); ++i) {
      payload[i] = payload[i] ^ mask[i % 4];
    }
  }

  return true;
}

void GHTTPConnection::ReceiverLoop() {
  // reads "frames" from the websocket. Currently 
  // there can be three different frames:
  // 0x1 - text frame
  // 0x2 - binary frame
  // 0x8 - close frame
  while(!fStopRequested) {
    unsigned char opcode = 0;
    std::vector<char> payload;
    if(!ReadFrame(opcode, payload)) {
      break;
    }

    if(opcode == 0x8) {
      break;
    }
    if(opcode == 0x1) {
      PrintTextFrame(payload);
    } else if(opcode == 0x2) {
      PrintSnapshotSummary(payload);
    }
  }

  fConnected = false;
}

void GHTTPConnection::PrintTextFrame(const std::vector<char>& payload) const {
  std::string text(payload.begin(), payload.end());
  std::istringstream stream(text);
  std::string line;
  while(std::getline(stream, line)) {
    if(!line.empty()) {
      std::cout << "\tlive text: " << line << std::endl;
    }
  }
}

void GHTTPConnection::PrintSnapshotSummary(const std::vector<char>& payload) const {
  const std::string separator = "\n\n";
  auto sep = std::search(payload.begin(), payload.end(),
                         separator.begin(), separator.end());
  if(sep == payload.end()) {
    std::cout << "\tlive binary snapshot: " << payload.size()
              << " bytes, no header found" << std::endl;
    return;
  }

  std::string header(payload.begin(), sep);
  std::string path;
  std::string class_name;
  std::string bytes;

  std::istringstream stream(header);
  std::string line;
  while(std::getline(stream, line)) {
    if(line.compare(0, 5, "PATH ") == 0) {
      path = line.substr(5);
    } else if(line.compare(0, 6, "CLASS ") == 0) {
      class_name = line.substr(6);
    } else if(line.compare(0, 6, "BYTES ") == 0) {
      bytes = line.substr(6);
    }
  }

  std::cout << "\tlive snapshot: " << path
            << " " << class_name
            << " " << bytes << " bytes" << std::endl;
}

void GHTTPConnection::CloseSocket() {
  int socket = fSocket;
  fSocket = -1;
  if(socket >= 0) {
    shutdown(socket, SHUT_RDWR);
    close(socket);
  }
}
