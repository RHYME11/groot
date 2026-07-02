#ifndef __GHTTPCONNECTION_H__
#define __GHTTPCONNECTION_H__

#include <atomic>
#include <string>
#include <thread>
#include <vector>

class GHTTPConnection {
public:
  GHTTPConnection();
  ~GHTTPConnection();

  bool Start(const std::string& url);
  void Stop();

  bool IsConnected() const { return fConnected; }

private:
  struct ParsedUrl {
    std::string host;
    std::string port;
    std::string path;
  };

  static bool ParseUrl(const std::string& url, ParsedUrl& parsed);
  static std::string Base64Encode(const unsigned char* data, size_t length);
  static std::string MakeWebSocketKey();

  bool ConnectSocket();
  bool SendHandshake();
  bool SendText(const std::string& text);
  bool ReadExact(void* data, size_t length);
  bool ReadFrame(unsigned char& opcode, std::vector<char>& payload);
  void ReceiverLoop();
  void PrintTextFrame(const std::vector<char>& payload) const;
  void PrintSnapshotSummary(const std::vector<char>& payload) const;
  void CloseSocket();

  std::string fUrl;
  ParsedUrl fParsedUrl;
  int fSocket;
  std::atomic<bool> fStopRequested;
  std::atomic<bool> fConnected;
  std::thread fThread;
};

#endif
