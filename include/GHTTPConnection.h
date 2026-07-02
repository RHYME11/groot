#ifndef __GHTTPCONNECTION_H__
#define __GHTTPCONNECTION_H__

#include <atomic>
#include <functional>
#include <string>
#include <thread>
#include <vector>

class TObject;

class GHTTPConnection {
public:
  using SnapshotCallback = std::function<void(const std::string& source,
                                              const std::string& path,
                                              TObject* object)>;

  GHTTPConnection();
  ~GHTTPConnection();

  void SetSnapshotCallback(SnapshotCallback callback);
  bool Start(const std::string& url);
  void Stop();

  bool IsConnected() const { return fConnected; }
  const std::string& SourceName() const { return fSourceName; }

private:
  struct ParsedUrl {
    std::string host;
    std::string port;
    std::string path;
  };

  struct Snapshot {
    std::string path;
    std::string className;
    std::vector<char> payload;
  };

  static bool ParseUrl(const std::string& url, ParsedUrl& parsed);
  static std::string Base64Encode(const unsigned char* data, size_t length);
  static std::string MakeWebSocketKey();
  static std::string MakeSourceName(const ParsedUrl& parsed);
  static std::string StripHistogramPrefix(const std::string& path);

  bool ConnectSocket();
  bool SendHandshake();
  bool SendText(const std::string& text);
  bool ReadExact(void* data, size_t length);
  bool ReadFrame(unsigned char& opcode, std::vector<char>& payload);
  bool ParseSnapshot(const std::vector<char>& frame, Snapshot& snapshot) const;
  TObject* DeserializeSnapshot(const Snapshot& snapshot) const;
  void ReceiverLoop();
  void PrintTextFrame(const std::vector<char>& payload) const;
  void PrintSnapshotSummary(const std::vector<char>& payload);
  void CloseSocket();

  std::string fUrl;
  ParsedUrl fParsedUrl;
  std::string fSourceName;
  int fSocket;
  std::atomic<bool> fStopRequested;
  std::atomic<bool> fConnected;
  std::thread fThread;
  SnapshotCallback fSnapshotCallback;
};

#endif
