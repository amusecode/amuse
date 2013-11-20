#!/bin/bash

function download_data {
    NAME=$1
    wget "http://www.cita.utoronto.ca/~iliev/rtwiki/lib/exe/fetch.php?id=tests5-7&cache=cache&media=rt:t7:${NAME}2.bin" -O "${NAME}2.bin"
    wget "http://www.cita.utoronto.ca/~iliev/rtwiki/lib/exe/fetch.php?id=tests5-7&cache=cache&media=rt:t7:${NAME}_add2.bin" -O "${NAME}_add2.bin"
}

download_data c2ray
download_data coral
download_data flash
download_data licorice
download_data susa
download_data zeus
