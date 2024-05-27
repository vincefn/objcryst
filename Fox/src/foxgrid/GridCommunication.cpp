#include "GridCommunication.h"

int generateUniqueID() {

    static bool initialized = false;
    static int currentID;
    
    if (!initialized) {
        std::srand(std::time(nullptr));
        currentID = std::rand() % 1000;

        initialized = true;
    }
    return currentID++;
}

GridCommunication::GridCommunication()
{
}
GridCommunication::~GridCommunication()
{

}
char GridCommunication::CalculateXORChecksum(const string& data) {
    char checksum = 0;
    for (char c : data) {
        checksum ^= c;
    }
    return checksum;
}
 vector<int> GridCommunication::getUsedPorts()
 {
    std::vector<int> localPorts;
    vector<string> result;

    //runs netstat
    try {
        string command = "netstat -an";        
        array<char, 1024> buffer; 
        unique_ptr<FILE, decltype(&PCLOSE)> pipe(POPEN(command.c_str(), "r"), PCLOSE);
        if (!pipe) {
            return localPorts;
        }
        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
            result.push_back(buffer.data());
        }        
    } catch (const exception& e) {
        cerr << "Exception: " << e.what() << std::endl;
        return localPorts;
    }
    
    /*
    Proto  Local Address          Foreign Address        State
    TCP    0.0.0.0:135            0.0.0.0:0              LISTENING
    TCP    0.0.0.0:445            0.0.0.0:0              LISTENING
    TCP    0.0.0.0:2853           0.0.0.0:0              LISTENING
    TCP    0.0.0.0:2854           0.0.0.0:0              LISTENING
    TCP    0.0.0.0:3389           0.0.0.0:0              LISTENING
    */

    //analyze netstat output and save ports
    for (const auto& line : result) {
        std::istringstream iss(line);
        std::string proto, localAddress, foreignAddress, state;

        if (!(iss >> proto >> localAddress >> foreignAddress >> state)) {
            continue; // skip header or malformed lines
        }

        size_t colonPos = localAddress.find_last_of(':');
        if (colonPos != std::string::npos) {
            std::string portStr = localAddress.substr(colonPos + 1);
            int port = std::stoi(portStr); // Convert port string to integer
            localPorts.push_back(port); // Store the port
        }
    }

    return localPorts;
}
long long GridCommunication::getTimeStampNanoSeconds() 
{
    auto now = std::chrono::high_resolution_clock::now();
    auto now_ms = std::chrono::time_point_cast<std::chrono::nanoseconds>(now).time_since_epoch().count();
    return now_ms;
}
unsigned int GridCommunication::getTimeStampMinutes()
{
    auto now = std::chrono::high_resolution_clock::now();
    auto now_m = std::chrono::time_point_cast<std::chrono::minutes>(now).time_since_epoch().count();
    return now_m;
}
long long GridCommunication::getTimeStampSeconds()
{
    auto now = std::chrono::high_resolution_clock::now();
    auto now_m = std::chrono::time_point_cast<std::chrono::seconds>(now).time_since_epoch().count();
    return now_m;
}
string GridCommunication::GenerateUniqueIdentifier(const string& userSpecificString) {
    
    auto now = std::chrono::high_resolution_clock::now();
    auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now).time_since_epoch().count();
    std::string dateTimeStr = std::to_string(now_ms);
    std::string combinedStr = dateTimeStr + userSpecificString;
        
    std::hash<std::string> hash_fn;
    size_t hash = hash_fn(combinedStr);  
    
    std::ostringstream oss;
    oss << std::setw(32) << std::setfill('0') << std::hex << hash; 
    std::string uniqueIdentifier = oss.str();
    
    if (uniqueIdentifier.length() > 32) {
        uniqueIdentifier = uniqueIdentifier.substr(0, 32);
    }

    return uniqueIdentifier;
}
bool GridCommunication::lastReadOK(wxSocketBase *socket, wxUint32 len)
{
    if(socket->LastReadCount()!=len) {
        WriteLogMessage("ERROR: last Read counts: "+to_string(socket->LastReadCount()) +" != " + to_string(len));
        return false;
    }
    if(socket->Error()==true) {
        WriteLogMessage("ERROR: last Read error: "+ to_string(socket->LastError()));
        return false;
    }
    return true;
}
bool GridCommunication::lastWriteOK(wxSocketBase *socket, wxUint32 len)
{
    if(socket->LastWriteCount()!=len) {
        WriteLogMessage("ERROR: last Write counts: "+to_string(socket->LastWriteCount()) +" != " + to_string(len));
        return false;
    }
    if(socket->Error()==true) {
        WriteLogMessage("ERROR: last Write error: "+ to_string(socket->LastError()));
        return false;
    }
    return true;
}
short GridCommunication::SendData(wxSocketBase *socket, long long msgID, const char* data, const wxUint32 dataLen)
{
    const wxUint32 maxChunkSize = 10*1024;
    
    //sending the header
    MessageHeader msgH;

    msgH.msgID = msgID;
    msgH.msgLen = dataLen;
    msgH.chunkLength = 0;
    msgH.chunkChecksum = 0;

    //socket.SaveState();

    WriteLogMessage("SendData start ");
    if(socket==0) {
        WriteLogMessage("ERROR: SendData: Socket=0");
        return 100;
    }
    if(!socket->IsConnected()) {
        WriteLogMessage("ERROR: SendData: Disconnected");
        return 101;
    }
    //socket->SetFlags(wxSOCKET_WAITALL|wxSOCKET_BLOCK);
    socket->SetFlags(wxSOCKET_WAITALL);
    
    socket->Write(&msgH, sizeof(MessageHeader));
    if(!lastWriteOK(socket, sizeof(MessageHeader))) {
        return 1;
    }

    wxUint32 remaining = dataLen;
    while (remaining > 0) {
        char message[maxChunkSize];   
        wxUint32 restChunkSize = std::min(remaining, static_cast<wxUint32>(maxChunkSize - sizeof(MessageHeader)));
        msgH.chunkLength = restChunkSize;

        memcpy(message, &msgH, sizeof(MessageHeader));
        memcpy(message + sizeof(MessageHeader), data, restChunkSize);

        //if (restChunkSize + sizeof(MessageHeader) < maxChunkSize) {
        //    memset(message + sizeof(MessageHeader) + restChunkSize, 0, maxChunkSize - (sizeof(MessageHeader) + maxChunkSize));
        //}

        socket->Write(message, restChunkSize + sizeof(MessageHeader));
        if(!lastWriteOK(socket, restChunkSize + sizeof(MessageHeader))) {
            return 0;
        }
        data += restChunkSize;
        remaining -= restChunkSize;
    }
    WriteLogMessage("SendData finished ");
    return 0;
}
short GridCommunication::ReadData(wxSocketBase *socket, long long &msgID, vector<char> &data)
{
    const wxUint32 maxChunkSize = 10*1024;
    MessageHeader msgH;

    socket->SetFlags(wxSOCKET_WAITALL);
    //socket->SetFlags(wxSOCKET_WAITALL|wxSOCKET_BLOCK);
    socket->Read(&msgH, sizeof(MessageHeader));

    if (!lastReadOK(socket, sizeof(MessageHeader))) {
        return 1;
    }
    if(msgH.chunkLength!=0) return 2; //in this stage we are expecting only the header!   
    msgID = msgH.msgID;

    data.clear();
    //data.reserve(msgH.msgLen);

    wxUint32 received = 0;
    while (received < msgH.msgLen) {
        char message[maxChunkSize];
        wxUint32 toReceive = std::min(static_cast<wxUint32>((msgH.msgLen - received) + sizeof(MessageHeader)), maxChunkSize);
        socket->Read(message, toReceive);
        if (!lastReadOK(socket, toReceive)) {
            return 3;
        }

        std::memcpy(&msgH, message, sizeof(MessageHeader));
        //long long t = message + sizeof(MessageHeader);


        data.insert(data.end(), message + sizeof(MessageHeader), message + sizeof(MessageHeader) + msgH.chunkLength);
        
        //std::memcpy(data + received, message + sizeof(MessageHeader), msgH.chunkLength);
        
        received += msgH.chunkLength;
    }   
    data.push_back('\0');
    return 0;
}
wxString GridCommunication::GetWorkingDir()
{
    return m_working_dir;
}
void GridCommunication::WriteLogMessage(wxString msg, wxString filename)
{
#ifdef ALLOW_GRID_LOGS
   wxString name;
#ifdef WIN32
   name = GetWorkingDir() + _T("\\") + filename;
#else
   name = GetWorkingDir() + _T("/") + filename;
#endif
   wxFile logfile(name, wxFile::write_append);
   if(logfile.IsOpened())
   {
      wxDateTime datetime = wxDateTime::Now();
      logfile.Write(datetime.Format(_T("%X ")) + msg + _T("\n"));
      logfile.Close();
   }
#endif
}
short GridCommunication::SendDataThread(wxSocketBase *socket, long long msgID, wxString text)
{
    socket->SetTimeout(10);
    socket->SetFlags(wxSOCKET_WAITALL | wxSOCKET_BLOCK);

    wxSocketOutputStream SocketOutputStream(*socket);
    //wxBufferedOutputStream BufferStream(SocketOutputStream);
    //BufferStream.Write(text.data(), text.length());
    //BufferStream.Sync();

    wxTextOutputStream toStream(SocketOutputStream);
    toStream.Write(text);

    return 0;
}
short GridCommunication::ReceiveDataThread(wxSocketBase *socket, long long &msgID, wxString &text)
{
    socket->SetTimeout(10);
    socket->SetFlags(wxSOCKET_WAITALL | wxSOCKET_BLOCK);

    wxStringOutputStream SOStream(&text);

    wxSocketInputStream SocketInputStream(*socket);
    wxTextInputStream tiStream(SocketInputStream);

    //SOStream.Write(tiStream.GetInputStream());
    //https://slideplayer.com/slide/14458020/

    return 0;
}
short GridCommunication::SendData2(wxSocketBase *socket, long long msgID, const char* data, wxUint32 const dataLen, bool thread)
{
    WriteLogMessage("SendData start");
    const wxUint32 maxChunkSize = 32*1024;
    MessageHeader msgH;

    msgH.msgID = msgID;
    msgH.msgLen = dataLen;
    msgH.chunkLength = 0;
    msgH.chunkChecksum = 0;

    if(thread) {
        socket->SetFlags(wxSOCKET_WAITALL|wxSOCKET_BLOCK);
    } else {
        socket->SetFlags(wxSOCKET_WAITALL);
    }

    socket->Write(&msgH, sizeof(MessageHeader));
    if(!lastWriteOK(socket, sizeof(MessageHeader))) {
        return 1;
    }

    wxUint32 remaining = dataLen;
    while (remaining > 0) {
        char message[maxChunkSize] = {};   
        wxUint32 restChunkSize = std::min(remaining, static_cast<wxUint32>(maxChunkSize - sizeof(MessageHeader)));
        msgH.chunkLength = restChunkSize;

        memcpy(message, &msgH, sizeof(MessageHeader));
        memcpy(message + sizeof(MessageHeader), data, restChunkSize);

        //if (restChunkSize + sizeof(MessageHeader) < maxChunkSize) {
        //    memset(message + sizeof(MessageHeader) + restChunkSize, 0, maxChunkSize - (sizeof(MessageHeader) + maxChunkSize));
        //}

        socket->Write(message, restChunkSize + sizeof(MessageHeader));
        if(!lastWriteOK(socket, restChunkSize + sizeof(MessageHeader))) {
            return 2;
        }
        data += restChunkSize;
        remaining -= restChunkSize;
        
        MessageHeader confirm;
        socket->Read(&confirm, sizeof(MessageHeader));
        
        if(!lastReadOK(socket, sizeof(MessageHeader))) {
            return 3;
        }
        if(confirm.msgID != msgID) {
            return 10;
        }

    }

    WriteLogMessage("SendData end");
    return 0;
    
}
short GridCommunication::ReceiveData2(wxSocketBase *socket, long long &msgID, vector<char> &data, bool thread)
{
    const wxUint32 maxChunkSize = 32*1024;
    MessageHeader msgH;

    //socket->SetFlags(wxSOCKET_WAITALL);

    if(thread) {
        socket->SetFlags(wxSOCKET_WAITALL|wxSOCKET_BLOCK);
    } else {
        socket->SetFlags(wxSOCKET_WAITALL);
    }

    socket->Read(&msgH, sizeof(MessageHeader));

    if (!lastReadOK(socket, sizeof(MessageHeader))) {
        return 1;
    }
    if(msgH.chunkLength!=0) return 2; //in this stage we are expecting only the header!   
    msgID = msgH.msgID;

    data.clear();
    data.reserve(msgH.msgLen);

    wxUint32 received = 0;
    while (received < msgH.msgLen) {
        char message[maxChunkSize] = {}; 
        wxUint32 toReceive = std::min(static_cast<wxUint32>((msgH.msgLen - received) + sizeof(MessageHeader)), maxChunkSize);
        socket->Read(message, toReceive);
        if (!lastReadOK(socket, toReceive)) {
            return 3;
        }

        std::memcpy(&msgH, message, sizeof(MessageHeader));
        data.insert(data.end(), message + sizeof(MessageHeader), message + (sizeof(MessageHeader) + msgH.chunkLength));
        received += msgH.chunkLength;

        socket->Write(&msgH, sizeof(MessageHeader));
        if(!lastWriteOK(socket, sizeof(MessageHeader))) {
            return 5;
        }
    }   
    data.push_back('\0');
    return 0;
}
