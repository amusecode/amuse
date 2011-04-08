#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>
#include <unistd.h>

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

jint capacity_out[10]; //size of header

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

MPI_Comm intercom;

void print_object(JNIEnv *env, const char *message, jobject object) {

	jclass object_class = (*env)->FindClass(env, "java/lang/Object");

	jmethodID toString_method = (*env)->GetMethodID(env, object_class,
			"toString", "()Ljava/lang/String;");

	jstring some_string = (jstring) (*env)->CallObjectMethod(env, object,
			toString_method);

	const char *str;

	str = (*env)->GetStringUTFChars(env, some_string, NULL);

	fprintf(stderr, "%s%s\n", message, str);

	(*env)->ReleaseStringUTFChars(env, some_string, str);
}

void update_buffers(JNIEnv *env) {

	jclass message_class = (*env)->FindClass(env, "ibis/amuse/AmuseMessage");

	jmethodID get_buffers_mid = (*env)->GetMethodID(env, message_class,
			"getByteBuffers", "(Z)[Ljava/nio/ByteBuffer;");

	jobjectArray byte_buffers = (jobjectArray) (*env)->CallObjectMethod(env,
			message_out, get_buffers_mid, JNI_FALSE);

	jobject header_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 0);
	header_out = (jint *) (*env)->GetDirectBufferAddress(env, header_bytes);

	jobject int_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 1);
	ints_out = (jint *) (*env)->GetDirectBufferAddress(env, int_bytes);
	capacity_out[HEADER_INT_COUNT_INDEX] = (*env)->GetDirectBufferCapacity(env,
			int_bytes) / SIZEOF_INT;

	jobject long_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 2);
	longs_out = (jlong *) (*env)->GetDirectBufferAddress(env, long_bytes);
	capacity_out[HEADER_LONG_COUNT_INDEX] = (*env)->GetDirectBufferCapacity(
			env, long_bytes) / SIZEOF_LONG;

	jobject float_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 3);
	floats_out = (jfloat *) (*env)->GetDirectBufferAddress(env, float_bytes);
	capacity_out[HEADER_FLOAT_COUNT_INDEX] = (*env)->GetDirectBufferCapacity(
			env, float_bytes) / SIZEOF_FLOAT;

	jobject double_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 4);
	doubles_out = (jdouble *) (*env)->GetDirectBufferAddress(env, double_bytes);
	capacity_out[HEADER_DOUBLE_COUNT_INDEX] = (*env)->GetDirectBufferCapacity(
			env, double_bytes) / SIZEOF_DOUBLE;

	jobject boolean_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 5);
	booleans_out = (jboolean *) (*env)->GetDirectBufferAddress(env,
			boolean_bytes);
	capacity_out[HEADER_BOOLEAN_COUNT_INDEX] = (*env)->GetDirectBufferCapacity(
			env, boolean_bytes) / SIZEOF_BOOLEAN;
}

void ensure_primitive_output_capacity(JNIEnv *env) {
	int i;

	jboolean ok = JNI_TRUE;

	for (i = HEADER_INT_COUNT_INDEX; i < HEADER_SIZE; i++) {
		if (header_out[i] > capacity_out[i]) {
			ok = JNI_FALSE;
		}
	}

	if (!ok) {
		fprintf(stderr, "increasing capacity\n");

		jclass message_class =
				(*env)->FindClass(env, "ibis/amuse/AmuseMessage");

		if (message_class == 0) {
			return;
		}

		jmethodID ensure_capacity_mid = (*env)->GetMethodID(env, message_class,
				"ensurePrimitiveCapacity", "()Z");

		if (ensure_capacity_mid == 0) {
			return;
		}

		//no need for result
		(*env)->CallBooleanMethod(env, message_out, ensure_capacity_mid);

		// call jni method ourselves (hope this works)
		update_buffers(env);
	}
}

JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_init(JNIEnv *env,
		jobject this_object, jstring code_name) {

	//needed?
	//dlopen("libmpi.dylib", RTLD_GLOBAL);

	char** empty;

	MPI_Init(0, &empty);

	MPI_Comm_spawn(
			"/Users/niels/workspace/amuse/src/amuse/community/bhtree/bhtree_worker",
			MPI_ARGV_NULL, 1, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercom,
			MPI_ERRCODES_IGNORE);

	fprintf(stderr, "mpi_ibis_worker.c: size of int = %d\n", sizeof(int));

}

/*
 * Class:     ibis_amuse_CommunityCode
 * Method:    setRequestMessage
 * Signature: (Libis/amuse/AmuseMessage;)V
 */JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_setRequestMessage(
		JNIEnv *env, jobject this_object, jobject request_message) {
	print_object(env, "setting request message: ", request_message);

	jclass message_class = (*env)->FindClass(env, "ibis/amuse/AmuseMessage");

	jmethodID get_buffers_mid = (*env)->GetMethodID(env, message_class,
			"getByteBuffers", "(Z)[Ljava/nio/ByteBuffer;");

	jobjectArray byte_buffers = (jobjectArray) (*env)->CallObjectMethod(env,
			request_message, get_buffers_mid, JNI_FALSE);

	jobject header_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 0);
	header_in = (jint *) (*env)->GetDirectBufferAddress(env, header_bytes);

	jobject int_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 1);
	ints_in = (jint *) (*env)->GetDirectBufferAddress(env, int_bytes);

	jobject long_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 2);
	longs_in = (jlong *) (*env)->GetDirectBufferAddress(env, long_bytes);

	jobject float_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 3);
	floats_in = (jfloat *) (*env)->GetDirectBufferAddress(env, float_bytes);

	jobject double_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 4);
	doubles_in = (jdouble *) (*env)->GetDirectBufferAddress(env, double_bytes);

	jobject boolean_bytes = (*env)->GetObjectArrayElement(env, byte_buffers, 5);
	booleans_in = (jboolean *) (*env)->GetDirectBufferAddress(env,
			boolean_bytes);
}

/*
 * Class:     ibis_amuse_CommunityCode
 * Method:    setResultMessage
 * Signature: (Libis/amuse/AmuseMessage;)V
 */JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_setResultMessage(
		JNIEnv *env, jobject this_object, jobject result_message) {

	print_object(env, "setting result message: ", result_message);

	message_out = (*env)->NewGlobalRef(env, result_message);

	update_buffers(env);
}

/*
 * Class:     ibis_amuse_CommunityCode
 * Method:    call
 * Signature: ()V
 */JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_call(JNIEnv *env,
		jobject this_object) {
	int result;
	int return_code;
	MPI_Status status;
	char error[100];
	jclass exception_class;
	int call_count = header_in[HEADER_CALL_COUNT_INDEX];
	int mpi_header[8];

	fprintf(stderr,
			"Amuse/Ibis/MPI native code: doing call, function ID = %d, count = %d\n",
			header_in[HEADER_FUNCTION_ID_INDEX], call_count);

	header_out[HEADER_CALL_ID_INDEX] = header_in[HEADER_CALL_ID_INDEX];
	header_out[HEADER_FUNCTION_ID_INDEX] = header_in[HEADER_FUNCTION_ID_INDEX];
	header_out[HEADER_CALL_COUNT_INDEX] = header_in[HEADER_CALL_COUNT_INDEX];

	//we do not support the "redirect output" function, simply pretend it worked...
	if (header_in[HEADER_FUNCTION_ID_INDEX] == 1141573512) {
		header_out[HEADER_INT_COUNT_INDEX] = 1;
		ensure_primitive_output_capacity(env);
		ints_out[0] = 0;
		return;
	}

	mpi_header[0] = header_in[HEADER_FUNCTION_ID_INDEX];
	mpi_header[1] = header_in[HEADER_CALL_COUNT_INDEX];
	mpi_header[2] = header_in[HEADER_DOUBLE_COUNT_INDEX] / call_count;
	mpi_header[3] = header_in[HEADER_INT_COUNT_INDEX] / call_count;
	mpi_header[4] = header_in[HEADER_FLOAT_COUNT_INDEX] / call_count;
	mpi_header[5] = header_in[HEADER_STRING_COUNT_INDEX] / call_count;
	mpi_header[6] = header_in[HEADER_BOOLEAN_COUNT_INDEX] / call_count;
	mpi_header[7] = header_in[HEADER_LONG_COUNT_INDEX] / call_count;

	fprintf(stderr, "sending header\n");
	MPI_Bcast(mpi_header, 8, MPI_INT, MPI_ROOT, intercom);


	if (header_in[HEADER_DOUBLE_COUNT_INDEX] > 0) {
		fprintf(stderr, "sending doubles\n");
		MPI_Bcast(doubles_in, header_in[HEADER_DOUBLE_COUNT_INDEX], MPI_DOUBLE,
				MPI_ROOT, intercom);
	}

	if (header_in[HEADER_INT_COUNT_INDEX] > 0) {
		fprintf(stderr, "sending ints\n");
		MPI_Bcast(ints_in, header_in[HEADER_INT_COUNT_INDEX], MPI_INT,
				MPI_ROOT, intercom);
	}

	if (header_in[HEADER_FLOAT_COUNT_INDEX] > 0) {
		fprintf(stderr, "sending floats\n");
		MPI_Bcast(floats_in, header_in[HEADER_FLOAT_COUNT_INDEX], MPI_FLOAT,
				MPI_ROOT, intercom);
	}

	if (header_in[HEADER_STRING_COUNT_INDEX] > 0) {
		//TODO: something with strings!
		MPI_Bcast(ints_in, 0, MPI_INT, MPI_ROOT, intercom);
		MPI_Bcast(booleans_in, 0, MPI_CHARACTER, MPI_ROOT, intercom);

		jclass exception_class = (*env)->FindClass(env, "java/lang/Exception");
		sprintf(error, "eep! a string!");
		if (exception_class != NULL) {
			(*env)->ThrowNew(env, exception_class, error);
		}
		return;
	}

	if (header_in[HEADER_BOOLEAN_COUNT_INDEX] > 0) {
		fprintf(stderr, "sending booleans");

		MPI_Bcast(booleans_in, header_in[HEADER_BOOLEAN_COUNT_INDEX],
				MPI_LOGICAL, MPI_ROOT, intercom);
	}

	if (header_in[HEADER_LONG_COUNT_INDEX] > 0) {
		fprintf(stderr, "sending longs");

		MPI_Bcast(longs_in, header_in[HEADER_LONG_COUNT_INDEX], MPI_INTEGER8,
				MPI_ROOT, intercom);
	}

	//now wait for result...

	fprintf(stderr, "Amuse/Ibis/MPI native code: waiting for result of function %d\n",
			header_in[HEADER_FUNCTION_ID_INDEX]);

	mpi_header[4] = 57575;

	result = MPI_Recv(mpi_header, 8, MPI_INTEGER4, 0, 999, intercom, &status);

	if (result != MPI_SUCCESS) {
		fprintf(stderr, "error in mpi call = %d\n", result);

	} else {
		fprintf(stderr, "received result %d %d %d %d %d %d %d %d\n", mpi_header[0], mpi_header[1], mpi_header[2], mpi_header[3], mpi_header[4], mpi_header[5], mpi_header[6], mpi_header[7]);
	}

	return_code = mpi_header[0]; //-1 and -2 for error

	header_out[HEADER_DOUBLE_COUNT_INDEX] = mpi_header[2] * call_count;
	header_out[HEADER_INT_COUNT_INDEX] = mpi_header[3] * call_count;
	header_out[HEADER_FLOAT_COUNT_INDEX] = mpi_header[4] * call_count;
	header_out[HEADER_STRING_COUNT_INDEX] = mpi_header[5] * call_count;
	header_out[HEADER_BOOLEAN_COUNT_INDEX] = mpi_header[6] * call_count;
	header_out[HEADER_LONG_COUNT_INDEX] = mpi_header[7] * call_count;

	print_object(env, "got header, ensuring capacity of ", message_out);

	//ensure_primitive_output_capacity(env);

	//print_object(env, "got header, receiving data for ", message_out);

//	if (header_out[HEADER_DOUBLE_COUNT_INDEX] > 0) {
//		fprintf(stderr, "Amuse/Ibis/MPI native code: receiving %d doubles\n",
//				header_out[HEADER_DOUBLE_COUNT_INDEX]);
//
//		result = MPI_Recv(&doubles_out, header_out[HEADER_DOUBLE_COUNT_INDEX],
//				MPI_DOUBLE, 0, 999, intercom, NULL);
//
//		if (result != MPI_SUCCESS) {
//			fprintf(stderr, "error in mpi call = %d\n", result);
//
//		}
//
//	}
//
//	fprintf(stderr, "size of MPI_INT = %d ", sizeof(MPI_INT) );
//
//	if (header_out[HEADER_INT_COUNT_INDEX] > 0) {
//		fprintf(stderr, "Amuse/Ibis/MPI native code: receiving %d ints\n",
//				header_out[HEADER_INT_COUNT_INDEX]);

		int array[1];

		result = MPI_Recv(array, 1, MPI_INTEGER4, 0, 999, intercom, &status);

		if (result != MPI_SUCCESS) {
					fprintf(stderr, "error in mpi call = %d\n", result);
				}

		fprintf(stderr, "Amuse/Ibis/MPI native code: received!! %d ints\n",
				header_out[HEADER_INT_COUNT_INDEX]);

//	}

	if (header_out[HEADER_FLOAT_COUNT_INDEX] > 0) {
		fprintf(stderr, "Amuse/Ibis/MPI native code: receiving %d floats\n",
				header_out[HEADER_FLOAT_COUNT_INDEX]);

		MPI_Recv(floats_out, header_out[HEADER_FLOAT_COUNT_INDEX], MPI_FLOAT,
				0, 999, intercom, NULL);
	}

	if (header_out[HEADER_STRING_COUNT_INDEX] > 0) {
		fprintf(stderr, "Amuse/Ibis/MPI native code: eep! strings!\n");

		//TODO: something with strings!
		return_code = -3;
	}

	if (header_out[HEADER_BOOLEAN_COUNT_INDEX] > 0) {
		fprintf(stderr, "Amuse/Ibis/MPI native code: receiving %d boolenas\n",
				header_out[HEADER_BOOLEAN_COUNT_INDEX]);

		MPI_Recv(booleans_out, header_out[HEADER_BOOLEAN_COUNT_INDEX],
				MPI_LOGICAL, 0, 999, intercom, NULL);
	}

	if (header_out[HEADER_LONG_COUNT_INDEX] > 0) {
		fprintf(stderr, "Amuse/Ibis/MPI native code: receiving %d longs\n",
				header_out[HEADER_LONG_COUNT_INDEX]);

		MPI_Recv(longs_out, header_out[HEADER_LONG_COUNT_INDEX], MPI_INTEGER8,
				0, 999, intercom, NULL);
	}

	fprintf(stderr, "checking for error in result\n");

	if (return_code != header_in[HEADER_FUNCTION_ID_INDEX]) {
		jclass exception_class = (*env)->FindClass(env, "java/lang/Exception");
		sprintf(error, "error %d in calling code function %id in call %d",
				return_code,
				header_in[HEADER_FUNCTION_ID_INDEX],
				header_in[HEADER_CALL_ID_INDEX]);

		if (exception_class != NULL) {
			(*env)->ThrowNew(env, exception_class, error);
		}
		return;
	}

	fprintf(stderr, "done!\n");

	///  comm.Recv(array,  source=0, tag=999)

	//	MPI_Comm_disconnect(&intercom);

}

