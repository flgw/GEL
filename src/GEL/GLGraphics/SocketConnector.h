//
// Created by flgw on 1/11/17.
//

#ifndef GEL_SOCKETCONNECTOR_H
#define GEL_SOCKETCONNECTOR_H

#include <iostream>
#include <sys/socket.h>
#include <sys/fcntl.h>
#include <unistd.h>
#include <cstring>
#include <mutex>
#include <thread>
#include <chrono>

namespace SocketConnector {
    std::string addr = "/tmp/MeshEdit";
    int sck = -1, sck_conn = -1;
    std::mutex mtx;
    std::string chainOfCommand;

    void open_socket() {

        // Create a socket with a local domain and a duplex stream type
        sck = socket(PF_LOCAL, SOCK_STREAM, 0);


        // Now, bind the socket to a point in the file system
        sockaddr sck_addr;
        sck_addr.sa_family = AF_LOCAL;
        memcpy(sck_addr.sa_data, addr.c_str(), addr.length());
#ifndef __GNUC__
        sck_addr.sa_len = addr.length();
#endif

        if (bind(sck, &sck_addr, sizeof(sockaddr)) != 0) {
            printf("Failed to bind socket");
            return;
        }

        // Listen for incoming connections (max 2)
        if (listen(sck, 2) != 0) {
            printf("Listening failed");
            return;
        }

        // Accept any takers
        sockaddr sck_addr_accept;
        socklen_t sck_len_accept;
        sck_conn = accept(sck, &sck_addr_accept, &sck_len_accept);
        if (sck_conn == -1) {
            printf("accept failed");
            return;
        }

        // set sockopts
        fcntl(sck_conn, F_SETFL, O_RDWR);

        // Send a greeting to the other end!
        std::string message = "MeshEdit socket connection\n\n";
        send(sck_conn, message.c_str(), message.length(), 0);
    }

    void close_socket() {
        close(sck);
    }

    bool listen_commands() {
        open_socket();

        char buffer[1024];
        std::string last_cmd;
        ssize_t l = 0;
        while ((l = recv(sck_conn, buffer, 1024, 0)) > 0) {

            /*std::cout << std::endl << "Size: " << l << std::endl << std::endl;
            auto str = std::string(buffer, l);
            std::cout << str.c_str() << std::endl;*/

            if (chainOfCommand.size() < 1024 * 1024) {
                mtx.lock();
                chainOfCommand.append(buffer, l);
                mtx.unlock();
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
        std::cout << "Shutting down socket" << std::endl;
        return false;
    }

    void execute_command(GLGraphics::Console *theConsole) {
        if (SocketConnector::mtx.try_lock()) {
            size_t pos1 = 0, pos2;

            for (pos2 = chainOfCommand.find('\n');
                 pos2 != std::string::npos;
                 pos1 = pos2 + 1, pos2 = chainOfCommand.find('\n', pos1)) {
                //std::cout << pos1 << " to " << pos2 << std::endl;
                std::cout << chainOfCommand.substr(pos1, pos2 - pos1) << std::endl;
                theConsole->execute(chainOfCommand.substr(pos1, pos2 - pos1).c_str());
            }

            chainOfCommand.erase(0, pos1);
            mtx.unlock();
        }
    }
};


#endif //GEL_SOCKETCONNECTOR_H
