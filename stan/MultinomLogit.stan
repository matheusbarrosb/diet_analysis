data {
int K;            // NUMBER OF CATEGORIES
int N;            // NUMBER OF OBSERVATIONS
int S;            // NUMBER OF SITES
int y[N];         // OBSERVATIONS
real TL[N];       // CONTINUOUS PREDICTOR. FISH SIZE
int site[N];      // CATEGORICAL PREDICTOR
}
parameters {
vector[K] b0[S];  // LOG-ODDS PREY PROBABILITIES
vector[K] bTL;    // EFFECT OF SIZE
}
transformed parameters {
matrix[N, K] LP;      // LINEAR PREDICTOR
for (i in 1:N) {
for (k in 1:K) {
for (s in 1:S) {
LP[i,k] = bTL[k]*TL[i] + b0[site[i],k];
}
}
}
}
model {
for (k in 1:K) {
for (s in 1:S) {
b0[s,k] ~ normal(0, 5);
}
}
bTL ~ normal(0, 1);
for (i in 1:N) {
y[i] ~ categorical_logit(LP[i]');
}
}
generated quantities {
vector[K] DC[S];
DC = rep_array(rep_vector(0, K), S);      // VECTOR FOR DIET COMPOSITION
for (i in 1:N) {
DC[site[i]] = softmax(bTL + b0[site[i]]);
}
for (k in 1:K) {
DC[k] /= sum(DC[k]);
}
}






