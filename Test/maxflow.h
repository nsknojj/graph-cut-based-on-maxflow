#ifndef MAXFLOW_H
#define MAXFLOW_H
#include <cstdio>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "maxflow.h"

void addedge(int u, int v, float c);

bool bfs();

float dfs(int u, float c);

void getflow();

void cut();

#endif