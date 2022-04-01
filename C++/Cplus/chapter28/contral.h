#pragma once
struct Squart {
	constexpr int operator()(int i){return i * i; }
};

struct Cube {
	constexpr int operator()(int i) { return i * i * i; }
};