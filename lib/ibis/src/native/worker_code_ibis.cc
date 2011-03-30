#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "bhtree_code.h"
#include "worker_code.h"
#include "stopcond.h"

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
const int SIZEOF_FLOAT = 4;
const int SIZEOF_DOUBLE = 8;
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

void update_buffers(JNIEnv *env) {

	jclass message_class = env->FindClass("ibis/amuse/AmuseMessage");

	jmethodID get_buffers_mid = env->GetMethodID(message_class,
			"getByteBuffers", "(Z)[Ljava/nio/ByteBuffer;");

	jobjectArray byte_buffers = (jobjectArray) env->CallObjectMethod(
			message_out, get_buffers_mid, false);

	jobject header_bytes = env->GetObjectArrayElement(byte_buffers, 0);
	header_out = (jint *) env->GetDirectBufferAddress(header_bytes);

	jobject int_bytes = env->GetObjectArrayElement(byte_buffers, 1);
	ints_out = (jint *) env->GetDirectBufferAddress(int_bytes);
	capacity_out[HEADER_INT_COUNT_INDEX] = env->GetDirectBufferCapacity(
			int_bytes) / SIZEOF_INT;

	jobject long_bytes = env->GetObjectArrayElement(byte_buffers, 2);
	longs_out = (jlong *) env->GetDirectBufferAddress(long_bytes);
	capacity_out[HEADER_LONG_COUNT_INDEX] = env->GetDirectBufferCapacity(
				long_bytes) / SIZEOF_LONG;

	jobject float_bytes = env->GetObjectArrayElement(byte_buffers, 3);
	floats_out = (jfloat *) env->GetDirectBufferAddress(float_bytes);
	capacity_out[HEADER_FLOAT_COUNT_INDEX] = env->GetDirectBufferCapacity(
				float_bytes) / SIZEOF_FLOAT;

	jobject double_bytes = env->GetObjectArrayElement(byte_buffers, 4);
	doubles_out = (jdouble *) env->GetDirectBufferAddress(double_bytes);
	capacity_out[HEADER_DOUBLE_COUNT_INDEX] = env->GetDirectBufferCapacity(
				double_bytes) / SIZEOF_DOUBLE;

	jobject boolean_bytes = env->GetObjectArrayElement(byte_buffers, 5);
	booleans_out = (jboolean *) env->GetDirectBufferAddress(boolean_bytes);
	capacity_out[HEADER_BOOLEAN_COUNT_INDEX] = env->GetDirectBufferCapacity(
				boolean_bytes) / SIZEOF_BOOLEAN;
}

void ensure_primitive_output_capacity(JNIEnv *env) {

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

		//no need for result
		env->CallBooleanMethod(message_out, ensure_capacity_mid);

		// call jni method ourselves (hope this works)
		update_buffers(env);
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

	jobject float_bytes = env->GetObjectArrayElement(byte_buffers, 3);
	floats_in = (jfloat *) env->GetDirectBufferAddress(float_bytes);

	jobject double_bytes = env->GetObjectArrayElement(byte_buffers, 4);
	doubles_in = (jdouble *) env->GetDirectBufferAddress(double_bytes);

	jobject boolean_bytes = env->GetObjectArrayElement(byte_buffers, 5);
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

	update_buffers(env);
}

/*
 * Class:     ibis_amuse_CommunityCode
 * Method:    call
 * Signature: ()V
 */JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_call(JNIEnv *env,
		jobject this_object) {
	int call_count = header_in[HEADER_CALL_COUNT_INDEX];
	fprintf(stderr, "doing call, function ID = %d\n",
			header_in[HEADER_FUNCTION_ID_INDEX]);

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
	case 1768994498:
		header_out[HEADER_INT_COUNT_INDEX] = 1;
		ensure_primitive_output_capacity(env);
		ints_out[0] = initialize_code();
		break;
	case 2069478464:
		header_out[HEADER_INT_COUNT_INDEX] = 1;
		ensure_primitive_output_capacity(env);
		ints_out[0] = commit_parameters();
		break;
	case 290264013:
		header_out[HEADER_INT_COUNT_INDEX] = (2 * call_count);
		ensure_primitive_output_capacity(env);

		for (int i = 0; i < call_count; i++) {
			ints_out[i] = new_particle(&ints_out[(1 * call_count) + i],
					doubles_in[i], doubles_in[(1 * call_count) + i],
					doubles_in[(2 * call_count) + i],
					doubles_in[(3 * call_count) + i],
					doubles_in[(4 * call_count) + i],
					doubles_in[(5 * call_count) + i],
					doubles_in[(6 * call_count) + i],
					doubles_in[(7 * call_count) + i]);
		}
		break;
	case 20920053:
		header_out[HEADER_INT_COUNT_INDEX] = 1;
		ensure_primitive_output_capacity(env);

		ints_out[0] = commit_particles();
		break;
	case 967950880:
		header_out[HEADER_INT_COUNT_INDEX] = (1 * call_count);
		header_out[HEADER_DOUBLE_COUNT_INDEX] = (8 * call_count);
		ensure_primitive_output_capacity(env);

		for (int i = 0; i < call_count; i++) {
			ints_out[i] = get_state(ints_in[i], &doubles_out[i],
					&doubles_out[(1 * call_count) + i],
					&doubles_out[(2 * call_count) + i],
					&doubles_out[(3 * call_count) + i],
					&doubles_out[(4 * call_count) + i],
					&doubles_out[(5 * call_count) + i],
					&doubles_out[(6 * call_count) + i],
					&doubles_out[(7 * call_count) + i]);
		}
		break;
	case 128926247:
		header_out[HEADER_INT_COUNT_INDEX] = 2;
		ensure_primitive_output_capacity(env);

		ints_out[0] = get_index_of_first_particle(&ints_out[1]);
		break;
	case 678380482:
		header_out[HEADER_INT_COUNT_INDEX] = 2;
		ensure_primitive_output_capacity(env);

		ints_out[0] = get_index_of_next_particle(ints_in[0], &ints_out[1]);
		break;
	case 1644113439:
		header_out[HEADER_INT_COUNT_INDEX] = 1;
		ensure_primitive_output_capacity(env);

		ints_out[0] = cleanup_code();
		break;
	default:
		jclass exception_class = env->FindClass(
				"java/lang/NoSuchMethodException");
		char error[100];
		sprintf(error, "unknown functon id %d in call %d",
				header_in[HEADER_FUNCTION_ID_INDEX],
				header_in[HEADER_CALL_ID_INDEX]);

		if (exception_class != NULL) {
			env->ThrowNew(exception_class, error);
		}
		return;
	}
}
