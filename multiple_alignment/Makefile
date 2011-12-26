CXX=g++
CXXFLAGS=-Wall -Werror -Wextra -I../overlap

multiple-alignment-test: multiple_alignment.o multiple_alignment_test.o multiple_alignment_builder.o ../overlap/overlapper.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm *.o 
