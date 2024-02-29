#include <curl/curl.h>
#include <string.h>
#include <stdio.h>

size_t write_callback(char *ptr, size_t size, size_t nmemb, void *userdata) {
    size_t total_size = size * nmemb;
    static size_t offset = 0;
    double *buffer = (double *)userdata;

    if (total_size > 0) {
        memcpy(buffer + offset, ptr, total_size);
        offset += total_size / sizeof(double);
    }

    return total_size;
}

int send_data_to_server(double* arr, int n) {
    CURL *curl = curl_easy_init();
    if (!curl) return -1;

    struct curl_slist *headers = curl_slist_append(NULL, "Content-Type: application/octet-stream");
    curl_easy_setopt(curl, CURLOPT_URL, "http://128.55.64.10:8080/receive_data");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, arr);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, sizeof(double) * n);
    CURLcode res = curl_easy_perform(curl);

    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
    return (res == CURLE_OK) ? 0 : -1;
}

int fetch_data_from_server(double* arr, int n) {
    CURL *curl = curl_easy_init();
    if (!curl) return -1;

    static size_t offset = 0;
    offset = 0;  // Reset offset for each fetch

    curl_easy_setopt(curl, CURLOPT_URL, "http://128.55.64.10:8080/send_data");
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, arr);
    CURLcode res = curl_easy_perform(curl);

    curl_easy_cleanup(curl);
    return (res == CURLE_OK) ? 0 : -1;
}
