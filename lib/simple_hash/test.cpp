#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern "C"
{
#include "simple_hash.h"
}

static const char* whitespace = " \t\r\n";

int main(int argc, const char* argv[])
{
    struct simple_hash ht;
    init_hash(&ht,128);
    for (;;)
    {
        char buf[256];
        if (fgets(buf, 256, stdin) == NULL)
            break;
        char *command = strtok(buf, whitespace);
        unsigned int key;
        unsigned int value;
        if (strcmp(command, "insert") == 0)
        {
            sscanf(strtok(NULL, whitespace), "%u", &key);
            sscanf(strtok(NULL, whitespace), "%u", &value);
            hash_insert(&ht,key,value);
        }
        else if (strcmp(command, "lookup") == 0)
        {
            sscanf(strtok(NULL, whitespace), "%u", &key);
            size_t value;
            int result;
            result=hash_lookup(&ht,key, &value);
            if (result==0)
                printf("%u\n", (unsigned int) value);
            else
                printf("None\n");
        }
        else if (strcmp(command, "increment") == 0)
        {
            sscanf(strtok(NULL, whitespace), "%u", &key);
            hash_cell_insert(&ht,key)->value++;
        }
        else if (strcmp(command, "delete") == 0)
        {
            sscanf(strtok(NULL, whitespace), "%u", &key);
            hash_delete(&ht,key);
        }
        else if (strcmp(command, "clear") == 0)
        {
            clear_hash(&ht);
        }
        else if (strcmp(command, "compact") == 0)
        {
            compact_hash(&ht);
        }
        fflush(stdout);
    }

    // Dump entire table
    printf("{\n");
    if(ht.m_zeroUsed) printf("    %u: %u,\n", ht.m_zeroCell.key, ht.m_zeroCell.value);
    for(size_t i=0;i<ht.m_arraySize;i++) if(ht.m_cells[i].key) printf("    %u: %u,\n", ht.m_cells[i].key, ht.m_cells[i].value);
    printf("}\n");
    end_hash(&ht);
    return 0;
}
