#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

const size_t BUF_SIZE = 1024;

struct thread_args_ {
    char* gp_out_path;

} thread_args_t;

void* dispatch(void* in) {

    thread_args_t* args = (thread_args_t*) in;

    int pipe[2];
    pipe(pipe);

    pid_t child = fork();

    // child process
    if (!child) {
        // execvp("gp.c", arg->gp_args);
        close(STDIN_FILENO);
        dup2(STDIN_FILENO, pipe[0]);
        system("gp -f -q -s 120000000 KH unpack_matrix.gp");
    }

    int fd = open(args->gp_out_path);
    char* buf[BUF_SIZE];
    
    ssize_t bytes_read;
    do {
        bytes_read = read(fd, buf, BUF_SIZE);
        
        ssize_t bytes_written = 0;
        while (bytes_written < bytes_read) {
            ssize_t ret = write(pipe[1], buf, bytes_read);
            if (ret < 0) {
                perror("write failed");
                break;
            }
            else {
                bytes_written += ret;
            }
        }        
    } while (bytes_read > 0);
    
    close(fd);
    close(pipe[1]);

    waitpid(child, 0, 0);
    return NULL;
}

int main() {

    size_t n_threads = 1;
    pthread_t tids[n_threads];

    for (size_t i = 0; i < n_threads; i++) {
        thread_args_t* t_args = malloc(sizeof(thread_args_t));
        asprintf(t_args->gp_out_path, "./gp_out%d", i);

        pthread_create(tids + i, NULL, dispatch, NULL);
    }

    for (size_t i = 0; i < n_threads; i++) {
        
        pthread_join(tids[i], NULL);
        
    }

    return 0;
}