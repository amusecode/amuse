#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "ibis_amuse_CommunityCode.h"


	/*
	 * Class:     ibis_amuse_CommunityCode
	 * Method:    setRequestMessage
	 * Signature: (Libis/amuse/AmuseMessage;)V
	 */
	JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_setRequestMessage(
			JNIEnv *env, jobject this_object, jobject request_message) {
		fprintf(stderr, "setting request message\n");

		jclass string_class = env->FindClass("java/lang/String");

		jclass this_class = env->FindClass("ibis/amuse/CommunityCode");

		jmethodID mid = env->GetMethodID(this_class, "toString", "()Ljava/lang/String;");

		jstring some_string = (jstring)env->CallObjectMethod(this_object, mid);

		const char *str;

		str = env->GetStringUTFChars(some_string, NULL);

		fprintf(stderr, "%s\n", str);

		return;
	}

	/*
	 * Class:     ibis_amuse_CommunityCode
	 * Method:    setResultMessage
	 * Signature: (Libis/amuse/AmuseMessage;)V
	 */
	JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_setResultMessage(
			JNIEnv *, jobject, jobject) {
		fprintf(stderr, "setting result message\n");
		return;
	}

	/*
	 * Class:     ibis_amuse_CommunityCode
	 * Method:    call
	 * Signature: ()V
	 */
	JNIEXPORT void JNICALL Java_ibis_amuse_CommunityCode_call(JNIEnv *, jobject) {
		fprintf(stderr, "doing call\n");
		return;
	}
