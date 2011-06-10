// public methods
void forsockets_init(int32_t port);
void forsockets_close();

void forsockets_receive_integers(int32_t *integers, int32_t length);
void forsockets_receive_longs(int64_t *longs, int32_t length);
void forsockets_receive_floats(float *floats, int32_t length);
void forsockets_receive_doubles(double *doubles, int32_t length);
void forsockets_receive_booleans(bool *booleans, int32_t length);
void forsockets_receive_string(char *string, int32_t length);

void forsockets_send_integers(int32_t *integers, int32_t length);
void forsockets_send_longs(int64_t *longs, int32_t length);
void forsockets_send_floats(float *floats, int32_t length);
void forsockets_send_doubles(double *doubles, int32_t length);
void forsockets_send_booleans(bool *booleans, int32_t length);
void forsockets_send_string(char *string, int32_t length);
