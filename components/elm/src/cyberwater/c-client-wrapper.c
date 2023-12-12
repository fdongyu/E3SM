// http_client_combined.c
#include <stdio.h>
#include <curl/curl.h>
#include <string.h>

// Callback function for writing received data
size_t write_callback(char* ptr, size_t size, size_t nmemb, void* userdata) {
    size_t total_size = size * nmemb;
    if (total_size <= 5 * sizeof(float)) {
        memcpy(userdata, ptr, total_size);
        return total_size;
    }
    return 0;
}

// Function to send data to server
//int send_data_to_server(float arr[], int n) {
int send_data_to_server(double arr[], int n) {
    CURL *curl;
    CURLcode res;
    struct curl_slist *headers = NULL;

    curl = curl_easy_init();
    if (!curl) {
        fprintf(stderr, "Failed to initialize curl\n");
        return -1;
    }

    headers = curl_slist_append(headers, "Content-Type: application/octet-stream");

    curl_easy_setopt(curl, CURLOPT_URL, "http://128.55.64.30:8080/send_data");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, arr);
    //curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, sizeof(float) * n);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, sizeof(double) * n);

    res = curl_easy_perform(curl);

    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK) {
        fprintf(stderr, "Failed to send data: %s\n", curl_easy_strerror(res));
        return -1;
    }

    return 0;
}

// Function to fetch data from server
int fetch_data_from_server(float arr[], int n) {
    CURL *curl;
    CURLcode res;

    curl = curl_easy_init();
    if (!curl) {
        fprintf(stderr, "Failed to initialize curl\n");
        return -1;
    }

    curl_easy_setopt(curl, CURLOPT_URL, "http://128.55.64.30:8080/get_data");
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, arr);

    res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK) {
        fprintf(stderr, "Failed to fetch data: %s\n", curl_easy_strerror(res));
        return -1;
    }

    return 0;
}
