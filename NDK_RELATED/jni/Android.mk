LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)
LOCAL_CFLAGS := -fvisibility=hidden -DDEBUG
LOCAL_MODULE    := fir_filter
LOCAL_SRC_FILES := ../../fir_filter.c
include $(BUILD_EXECUTABLE)