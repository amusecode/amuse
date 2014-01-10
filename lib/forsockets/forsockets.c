#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#ifdef WIN32
	#include <winsock2.h>
#else
	#include <sys/socket.h>
	#include <netinet/in.h>
	#include <netdb.h>
	#include <netinet/tcp.h>
#endif

int32_t socketfd;

//private funtions

void forsockets_send(void *buffer, int32_t length, int32_t file_descriptor) {
	int32_t total_written = 0;
	int32_t written;

	while (total_written < length) {
		
#ifdef WIN32
		written = send(file_descriptor, ((char *) buffer) + total_written,
                        length - total_written, 0);
#else
		written = write(file_descriptor, ((char *) buffer) + total_written,
				length - total_written);
#endif

		if (written == -1) {
			fprintf(stderr, "could not write data\n");
			exit(1);
		}

		total_written = total_written + written;
	}

}

void forsockets_receive(void *buffer, int32_t length, int32_t file_descriptor) {
	int32_t total_read = 0;
	int32_t bytes_read;

	while (total_read < length) {
#ifdef WIN32
		bytes_read = recv(file_descriptor, ((char *) buffer) + total_read,
				length - total_read, 0);
#else
		bytes_read = read(file_descriptor, ((char *) buffer) + total_read,
				length - total_read);
#endif

//		fprintf(stderr, "received %d bytes, " + bytes_read);

		if (bytes_read == -1) {
			fprintf(stderr, "could not read data\n");
			exit(1);
		}

		total_read = total_read + bytes_read;
	}

}

//public functions

void forsockets_receive_integers(int32_t *integers, int32_t length) {
	forsockets_receive((void *) integers, length * sizeof(int32_t), socketfd);
}

void forsockets_receive_longs(int64_t *longs, int32_t length) {
	forsockets_receive((void *) longs, length * sizeof(int64_t), socketfd);
}

void forsockets_receive_floats(float *floats, int32_t length) {
	forsockets_receive((void *) floats, length * sizeof(float), socketfd);
}

void forsockets_receive_doubles(double *doubles, int32_t length) {
//	int i;
//	double array[100];
//	forsockets_receive((void *) array, length * sizeof(double), socketfd);
//
//	fprintf(stderr, "received doubles: ");
//	for (i = 0; i < length; i++) {
//		fprintf(stderr, " %e", array[i]);
//		doubles[i] = array[i];
//	}
//	fprintf(stderr, "\n");
	forsockets_receive((void *) doubles, length * sizeof(double), socketfd);
}

void forsockets_receive_booleans(bool *booleans, int32_t length) {
	forsockets_receive((void *) booleans, length * sizeof(bool), socketfd);
}

void forsockets_receive_string(char *string, int32_t length) {
	forsockets_receive((void *) string, length * sizeof(char), socketfd);
}

void forsockets_send_integers(int32_t *integers, int32_t length) {
	forsockets_send((void *) integers, length * sizeof(int32_t), socketfd);
}

void forsockets_send_longs(int64_t *longs, int32_t length) {
	forsockets_send((void *) longs, length * sizeof(int64_t), socketfd);
}

void forsockets_send_floats(float *floats, int32_t length) {
	forsockets_send((void *) floats, length * sizeof(float), socketfd);
}

void forsockets_send_doubles(double *doubles, int32_t length) {
	forsockets_send((void *) doubles, length * sizeof(double), socketfd);
}

void forsockets_send_booleans(bool *booleans, int32_t length) {
	forsockets_send((void *) booleans, length * sizeof(bool), socketfd);
}

void forsockets_send_string(char *string, int32_t length) {
	forsockets_send((void *) string, length * sizeof(char), socketfd);
}

void forsockets_init(int32_t port) {
	struct sockaddr_in serv_addr;
	struct hostent *server;
	int on = 1;
	
#ifdef WIN32
	WSADATA wsaData;
	int iResult;

	// Initialize Winsock
	iResult = WSAStartup(MAKEWORD(2,2), &wsaData);
	if (iResult != 0) {
	printf("WSAStartup failed: %d\n", iResult);
	exit(1);
	}
#endif

//	fprintf(stderr, "initializing forsockets\n");
//
//	fprintf(stderr, "sizeof float = %u\n", sizeof(float));
//	fprintf(stderr, "sizeof double = %u\n", sizeof(double));

	socketfd = socket(AF_INET, SOCK_STREAM, 0);

	if (socketfd < 0) {
		fprintf(stderr, "cannot open socket\n");
		exit(0);
	}

	server = gethostbyname("127.0.0.1");

//	fprintf(stderr, "connecting...\n");

	memset((char *) &serv_addr, '\0', sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	memcpy((char *) &serv_addr.sin_addr.s_addr, (char *) server->h_addr, 
			server->h_length);
	serv_addr.sin_port = htons(port);
	if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr))
			< 0) {
		fprintf(stderr, "cannot connect socket\n");
		exit(0);

	}
	
        setsockopt(socketfd, IPPROTO_TCP, TCP_NODELAY, (const char *) &on, sizeof(on));    

	/*fprintf(stderr, "finished initializing forsockets\n");*/
}

void forsockets_close() {
#ifdef WIN32
	closesocket(socketfd);
#else
	close(socketfd);
#endif
}
