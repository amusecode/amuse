#!/bin/bash

function download_data {
    NAME=$1
    wget "http://www.cita.utoronto.ca/~iliev/rtwiki/lib/exe/fetch.php?id=tests5-7&cache=cache&media=rt:t5:${NAME}1.bin" -O "${NAME}1.bin"
    wget "http://www.cita.utoronto.ca/~iliev/rtwiki/lib/exe/fetch.php?id=tests5-7&cache=cache&media=rt:t5:${NAME}4.bin" -O "${NAME}4.bin"
    wget "http://www.cita.utoronto.ca/~iliev/rtwiki/lib/exe/fetch.php?id=tests5-7&cache=cache&media=rt:t5:${NAME}5.bin" -O "${NAME}5.bin"
    wget "http://www.cita.utoronto.ca/~iliev/rtwiki/lib/exe/fetch.php?id=tests5-7&cache=cache&media=rt:t5:${NAME}_add1.bin" -O "${NAME}_add1.bin"
    wget "http://www.cita.utoronto.ca/~iliev/rtwiki/lib/exe/fetch.php?id=tests5-7&cache=cache&media=rt:t5:${NAME}_add4.bin" -O "${NAME}_add4.bin"
    wget "http://www.cita.utoronto.ca/~iliev/rtwiki/lib/exe/fetch.php?id=tests5-7&cache=cache&media=rt:t5:${NAME}_add5.bin" -O "${NAME}_add5.bin"
}

download_data susa
download_data rh1d
download_data zeus