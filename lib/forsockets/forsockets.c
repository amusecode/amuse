#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

int32_t socketfd;

//private funtions

void forsockets_send(void *buffer, int32_t length, int32_t file_descriptor) {
	int32_t total_written = 0;
	int32_t written;

	while (total_written < length) {
		written = write(file_descriptor, ((char *) buffer) + total_written,
				length - total_written);

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
		bytes_read = read(file_descriptor, ((char *) buffer) + total_read,
				length - total_read);

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

//	fprintf(stderr, "initializing forsockets\n");
//
//	fprintf(stderr, "sizeof float = %u\n", sizeof(float));
//	fprintf(stderr, "sizeof double = %u\n", sizeof(double));

	socketfd = socket(AF_INET, SOCK_STREAM, 0);

	if (socketfd < 0) {
		fprintf(stderr, "cannot open socket\n");
		exit(0);
	}

	server = gethostbyname("localhost");

//	fprintf(stderr, "connecting...\n");

	bzero((char *) &serv_addr, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr,
			server->h_length);
	serv_addr.sin_port = htons(port);
	if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr))
			< 0) {
		fprintf(stderr, "cannot connect socket\n");
		exit(0);

	}

	fprintf(stderr, "finished initializing forsockets\n");
}

void forsockets_close() {
	close(socketfd);
}
