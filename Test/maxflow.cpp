#include "maxflow.h"
using namespace std;

const int f_maxn = 2250010,
f_maxm = 13000000,
maxnum = 2147483647;
struct edge{
	int p, n, v;
	float c;
}	e[f_maxm];
int S, T, tot = 2;
int fir[f_maxn], he[f_maxn];
int list[f_maxn];
bool f[f_maxn];

void addedge(int u, int v, float c) {
	e[tot].p = v;
	e[tot].c = c;
	e[tot].n = fir[u];
	fir[u] = tot++;
}

bool bfs() {
	int h, r, k, u;
	memset(he, 0, sizeof(he));
	list[0] = T; he[T] = 1;
	for (h = r = 0; h <= r; h++){
		u = list[h];
		for (k = fir[u]; k; k = e[k].n)
		if ((e[k ^ 1].c > 0) && (!he[e[k].p])){
			list[++r] = e[k].p;
			he[e[k].p] = he[u] + 1;
		}
	}
	return he[S] != 0;
}

float dfs(int u, float c) {
	if (u == T) return c;
	float change, re;
	int k;
	for (re = 0, k = fir[u]; k; k = e[k].n) {
		if (e[k].c&&he[e[k].p] == he[u] - 1) {
			change = dfs(e[k].p, min(c, e[k].c));
			e[k].c -= change;
			e[k ^ 1].c += change;
			c -= change;
			re += change;
		}
		if (!c) break;
	}
	if (!re)he[u] = maxnum;
	return re;
}

void go(int x) {
	f[x] = 1;
	int l, r;
	l = r = 1;
	list[r] = x;
	for (; l <= r; l++) {
		x = list[l];
		for (int k = fir[x]; k; k = e[k].n)
		if (e[k].c > 0 && !f[e[k].p]){
			f[e[k].p] = 1;
			list[++r] = e[k].p;
		}
	}
}

void cut() {
	memset(f, 0, sizeof f);
	go(S);
}

void getflow() {
	while ( bfs() )
		if (dfs(S, maxnum) < 0.0000001) {};// break;
}
