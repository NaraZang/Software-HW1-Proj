#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define INFINITY 1e9 // Define a large value for infinity

/**
  Declaring Functions :
 */  

float** initialize_centroids(float** points, int K, int d, int N);
float** load_points(const char* filename, int* N_out, int* d_out);
float dist(float* p , float* q, int d);
void print_centroids(float **Centroids, int K, int d);
int* assign_clusters(float** points, float** centroids, int K, int d, int N);
int has_converged(float** old_centroids, float** new_centroids, int K, int d, float epsilon);
float** update_centroids(float** points, int* cluster_indices, int num_points, int K, int dim);
float** create_zero_centroids(int k, int d); //helper function Not used directly in K means



float** load_points(const char* filename, int* N_out, int* d_out) { //updated Temporarily reading from file
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        perror("File open failed");
        exit(1);
    }

    char line[MAX_LINE_LEN];   // Buffer for reading each line
    int capacity = 10;         // Initial capacity for number of points
    int N = 0;                 // Actual number of points read
    int d = -1;                // Number of dimensions (to be determined)

    // Allocate initial array of point pointers
    float** points = malloc(capacity * sizeof(float*));
    if (!points) {
        perror("malloc failed");
        exit(1);
    }

    // Read each line from the file
    while (fgets(line, MAX_LINE_LEN, fp)) {
        // Resize array if needed
        if (N >= capacity) {
            capacity *= 2;
            points = realloc(points, capacity * sizeof(float*));
            if (!points) {
                perror("realloc failed");
                exit(1);
            }
        }

        // Count how many numbers (dimensions) are in the line
        int curr_d = 0;
        char* tmp = strdup(line);  // Make a copy of the line to tokenize safely
        char* token = strtok(tmp, " \t\n");
        while (token) {
            curr_d++;
            token = strtok(NULL, " \t\n");
        }
        free(tmp);

        // Check consistency of dimension count
        if (d == -1) {
            d = curr_d;  // First line sets the dimensionality
        } else if (curr_d != d) {
            fprintf(stderr, "Inconsistent number of dimensions in line %d\n", N + 1);
            exit(1);
        }

        // Allocate memory for this point
        points[N] = malloc(d * sizeof(float));
        if (!points[N]) {
            perror("malloc failed");
            exit(1);
        }

        // Parse and store the float values
        int i = 0;
        token = strtok(line, " \t\n");
        while (token && i < d) {
            points[N][i++] = atof(token);
            token = strtok(NULL, " \t\n");
        }

        N++;  // One more point loaded
    }

    fclose(fp);  // Close the file

    // Set output values
    *N_out = N;
    *d_out = d;

    // Example: print the points - could delete later
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < d; j++) {
            printf("%.4f ", points[i][j]);
        }
        printf("\n");
    }
    return points; // Return the array 

    // Free memory
    /*
    The function that allocates the memory is not necessarily the one that frees it — 
    but someone must take responsibility to free it once and only once, when it’s no longer needed.
    
    for (int i = 0; i < K; i++) {
        free(points[i]);
    }
    free(points);
*/ 
}
float** initialize_centroids(float** points, int K, int d, int N) {
    float **copy = malloc(K * sizeof(float*));

    for (int i = 0; i < K; i++) {
        copy[i] = malloc(d * sizeof(float));
        for (int j = 0; j < d; j++) {
            copy[i][j] = points[i][j];
        }
    }

    return copy;
}
float dist(float* p, float* q , int d) {
    float sum = 0.0;
    for (int i = 0; i < d; i++) { 
        float diff = p[i] - q[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}
void print_centroids(float **Centroids, int K, int d){
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < d; j++) {
            printf("%.4f", Centroids[i][j]);
            if (j < d - 1)
                printf(","); // Add comma between values
        }
        printf("\n");
    }
}
int* assign_clusters(float** points, float** centroids, int K, int d, int N){
    int* cluster_indices = malloc(N * sizeof(int));  // ← correct size
    if (!cluster_indices) { // If malloc fails and returns NULL, we handle it and exit
        perror("malloc failed");
        exit(1);
    }
    for (int i=0 ; i< N ; i++){
        float min_dist = INFINITY;
        int min_index = -1; // check if better to name closest_centroid 
        for(int j = 0; j < K; j++) {
            float distance = dist(points[i], centroids[j], d);
            if (distance < min_dist) {
                min_dist = distance;
                min_index = j;
            }
        }
        cluster_indices[i] = min_index; // Assign the closest centroid index to the point
    }
    return cluster_indices;
}

float** update_centroids(float** points, int* cluster_indices, int K, int d, int N) {
    float** new_centroids=create_zero_centroids(K,d); //creating an empty (0.0) list of lists that each inner list of size dim
    int* counts = calloc(K, sizeof(int)); //allocating memory to count the number of points for each centroid
    for (i=0; i<N; i++) { 
        int cluster = cluster_indices[i]; // cluster is the index of the centroid to which the point belongs
        for (int s = 0; s < d; s++) {
            new_centroids[cluster][s] += points[i][s];} // Add the point's coordinates to the corresponding centroid - calculating the new centroid
        counts[cluster] += 1; // Increment the count for this centroid
    }
    for (int j = 0; j < K; j++) {
        if (counts[j] == 0) {
            continue;  // Avoid division by zero 
        }
        for (int s = 0; s < d; s++) {
            new_centroids[j][s] /= counts[j];  // Average each dimension
            }
    }
    return new_centroids; // Return the new centroids
}


//create_zero_centroides is a helper function,used to create an empty (0.0) list of lists that each inner list of size d, namely new centroids

float** create_zero_centroids(int k, int d) { 
    float** centroids = malloc(k * sizeof(float*)); //allocates memory for K pointers to float,These pointers will each point to the start of a row.
    if (!centroids) {
        perror("malloc failed for centroids");
        return NULL;} 
    float* data = calloc(k * d, sizeof(float)); // zero-initializes all floats, uses calloc to allocate and initialize all bytes to zero.
    if (!data) {
        perror("calloc failed for data");
        free(centroids);
        return NULL; }
    for (int i = 0; i <k; i++) {
    centroids[i] = data + i * d; } //each elem in centroids is a pointer to a 0.0 float in data,Sets centroids[i] (the pointer to row i) to point to the ith chunk of dim floats inside the big block.
    return centroids; }

int has_converged(float** old_centroids, float** new_centroids, int K, int d, float epsilon){
    for (int i = 0; i < K; i++) {
        float distance = dist(old_centroids[i], new_centroids[i], d);
        if (distance >= epsilon) { //check if it should be greater or (greater or equal)
            return 0; // Not converged
        }
    }
    return 1; // Converged
    // in this case 1 returns true and 0 returns false
}
void k_means(const char* filename, int K, int iter) { //temporarily reading ftom file
    double epsilon = 0.001;
    int N, dim;

    float **points = load_points(filename,&N, &dim); //update load points
    if (!(1 < K && K < N)) {
        fprintf(stderr, "Invalid number of clusters!\n");
        exit(1);
    }
    if (!(1 < iter && iter < 1000)) {
        fprintf(stderr, "Invalid maximum iteration!\n");
        exit(1);
    }

    float **centroids = initialize_centroids(points, K, N, dim);
    int *cluster_indices;
    float **new_centroids;
    int i;

    for (i = 0; i < iter; i++) {
        cluster_indices = assign_clusters(points, centroids, K, dim, N);
        new_centroids = update_centroids(points, cluster_indices, K, dim, N);

        if (has_converged(centroids, new_centroids, K, dim, epsilon)) {
            break;
        }

        centroids = new_centroids;

        print_centroids(centroids, K, dim);
    }
    // TODO: free allocated memory here
}
void k_means_default(int K) { //MAYBE CORRECTED!!
    k_means(K, 400);  // call with default iter = 400
}



        
        




