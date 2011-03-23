#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "ibis_amuse_CommunityCode.h"

const int HEADER_SIZE = 10; // integers

const int HEADER_CALL_ID_INDEX = 1;
const int HEADER_FUNCTION_ID_INDEX = 2;
const int HEADER_CALL_COUNT_INDEX = 3;
const int HEADER_INT_COUNT_INDEX = 4;
const int HEADER_LONG_COUNT_INDEX = 5;
const int HEADER_FLOAT_COUNT_INDEX = 6;
const int HEADER_DOUBLE_COUNT_INDEX = 7;
const int HEADER_BOOLEAN_COUNT_INDEX = 8;
const int HEADER_STRING_COUNT_INDEX = 9;

const int SIZEOF_INT = 4;
const int SIZEOF_LONG = 8;
const int SIZEOF_FLOAT = 8;
const int SIZEOF_DOUBLE = 4;
const int SIZEOF_BOOLEAN = 1;

jint *header_in;

jint *ints_in;
jlong *longs_in;
jfloat *floats_in;
jdouble *doubles_in;
jboolean *booleans_in;
//TODO: something with strings/chars

jobject message_out;

jint *header_out;

jint *ints_out;
jlong *longs_out;
jfloat *floats_out;
jdouble *doubles_out;
jboolean *booleans_out;

jint *capacity_out = new jint[HEADER_SIZE];

void print_object(JNIEnv *env, const char *message, jobject object) {

	jclass object_class = env->FindClass("java/lang/Object");

	jmethodID toString_method = env->GetMethodID(object_class, "toString",
			"()Ljava/lang/String;");

	jstring some_string = (jstring) env->CallObjectMethod(object,
			toString_method);

	const char *str;

	str = env->GetStringUTFChars(some_string, NULL);

	fprintf(stderr, "%s%s\n", message, str);

	env->ReleaseStringUTFChars(some_string, str);
}

void ensure_primitive_output_capacity(JNIEnv *env) {

	print_object(env, "output message saved ", message_out);

	jboolean ok = true;

	for (int i = HEADER_INT_COUNT_INDEX; i < HEADER_SIZE; i++) {
		if (header_out[i] > capacity_out[i]) {
			ok = false;
		}
	}

	if (!ok) {
		fprintf(stderr, "increasing capacity\n");

		jclass message_class = env->FindClass("ibis/amuse/AmuseMessage");

		if (message_class == 0) {
			return;
		}

		jmethodID ensure_capacity_mid = env->GetMethodID(message_class,
				"ensurePrimitiveCapacity", "()Z");

		if (ensure_capacity_mid == 0) {
			return;
		}

		jboolean result = env->CallBooleanMethod(message_out, ensure_capacity_mid);

		// call jni method ourselves (hope this works)
		Java_ibis_amuse_CommunityCode_setResultMessage(env, message_out,
				message_out);
	}

}

/*
 * Class:     ibis_amuse_CommunityCode
 * Method:    setRequestMessage
 * Signature: (Libis/amuse/AmuseMessage;)V
 */JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_setRequestMessage(
		JNIEnv *env, jobject this_object, jobject request_message) {
	print_object(env, "setting request message: ", request_message);

	jclass message_class = env->FindClass("ibis/amuse/AmuseMessage");

	jclass bytebuffer_array_class = env->FindClass("[Ljava/nio/ByteBuffer;");

	jmethodID get_buffers_mid = env->GetMethodID(message_class,
			"getByteBuffers", "(Z)[Ljava/nio/ByteBuffer;");

	jobjectArray byte_buffers = (jobjectArray) env->CallObjectMethod(
			request_message, get_buffers_mid, false);

	jobject header_bytes = env->GetObjectArrayElement(byte_buffers, 0);
	header_in = (jint *) env->GetDirectBufferAddress(header_bytes);

	jobject int_bytes = env->GetObjectArrayElement(byte_buffers, 1);
	ints_in = (jint *) env->GetDirectBufferAddress(int_bytes);

	jobject long_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	longs_in = (jlong *) env->GetDirectBufferAddress(long_bytes);

	jobject float_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	floats_in = (jfloat *) env->GetDirectBufferAddress(float_bytes);

	jobject double_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	doubles_in = (jdouble *) env->GetDirectBufferAddress(double_bytes);

	jobject boolean_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	booleans_in = (jboolean *) env->GetDirectBufferAddress(boolean_bytes);
}

/*
 * Class:     ibis_amuse_CommunityCode
 * Method:    setResultMessage
 * Signature: (Libis/amuse/AmuseMessage;)V
 */JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_setResultMessage(
		JNIEnv *env, jobject this_object, jobject result_message) {
	print_object(env, "setting result message: ", result_message);

	message_out = env->NewGlobalRef(result_message);

	jclass message_class = env->FindClass("ibis/amuse/AmuseMessage");

	jclass bytebuffer_array_class = env->FindClass("[Ljava/nio/ByteBuffer;");

	jmethodID get_buffers_mid = env->GetMethodID(message_class,
			"getByteBuffers", "(Z)[Ljava/nio/ByteBuffer;");

	jobjectArray byte_buffers = (jobjectArray) env->CallObjectMethod(
			result_message, get_buffers_mid, false);

	jobject header_bytes = env->GetObjectArrayElement(byte_buffers, 0);
	header_out = (jint *) env->GetDirectBufferAddress(header_bytes);


	jobject int_bytes = env->GetObjectArrayElement(byte_buffers, 1);
	ints_out = (jint *) env->GetDirectBufferAddress(int_bytes);
	capacity_out[HEADER_INT_COUNT_INDEX] = env->GetDirectBufferCapacity(
			int_bytes) / SIZEOF_INT;

	jobject long_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	longs_out = (jlong *) env->GetDirectBufferAddress(long_bytes);

	jobject float_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	floats_out = (jfloat *) env->GetDirectBufferAddress(float_bytes);

	jobject double_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	doubles_out = (jdouble *) env->GetDirectBufferAddress(double_bytes);

	jobject boolean_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	booleans_out = (jboolean *) env->GetDirectBufferAddress(boolean_bytes);
}

/*
 * Class:     ibis_amuse_CommunityCode
 * Method:    call
 * Signature: ()V
 */JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_call(JNIEnv *env, jobject this_object) {
	fprintf(stderr, "doing call, function ID = %d\n", header_in[HEADER_FUNCTION_ID_INDEX]);

	header_out[HEADER_CALL_ID_INDEX] = header_in[HEADER_CALL_ID_INDEX];
	header_out[HEADER_FUNCTION_ID_INDEX] = header_in[HEADER_FUNCTION_ID_INDEX];
	header_out[HEADER_CALL_COUNT_INDEX] = header_in[HEADER_CALL_COUNT_INDEX];

	switch (header_in[HEADER_FUNCTION_ID_INDEX]) {
	case 0:
		//IGNORE
		break;
	case 1141573512:
		header_out[HEADER_INT_COUNT_INDEX] = 1;
		ensure_primitive_output_capacity(env);
		ints_out[0] = 0;
		break;
	}

	return;
}
