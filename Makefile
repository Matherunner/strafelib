CXX ?= g++
CXXFLAGS = -std=c++14 -flto -Wall -Wextra -Ofast -march=native -mtune=native
OUTPUT = test_strafelib
TEST_OBJS = test_strafelib.o

test: $(OUTPUT)

$(OUTPUT): $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(TEST_OBJS) $(OUTPUT)
