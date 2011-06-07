// public methods
int32_t forsockets_init(int32_t port);
int32_t forsockets_close();

int32_t forsockets_receive_ints(int32_t *ints, int32_t length);
int32_t forsockets_receive_longs(int64_t *longs, int32_t length);
int32_t forsockets_receive_floats(float *floats, int32_t length);
int32_t forsockets_receive_doubles(double *doubles, int32_t length);
int32_t forsockets_receive_booleans(bool *booleans, int32_t length);
int32_t forsockets_receive_string(char *string, int32_t length);

int32_t forsockets_send_ints(int32_t *ints, int32_t length);
int32_t forsockets_send_longs(int64_t *longs, int32_t length);
int32_t forsockets_send_floats(float *floats, int32_t length);
int32_t forsockets_send_doubles(double *doubles, int32_t length);
int32_t forsockets_send_booleans(bool *booleans, int32_t length);
int32_t forsockets_send_string(char *string, int32_t length);
