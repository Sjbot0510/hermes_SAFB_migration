# SAFB C Project — Build System
#
# Uses the $CC environment variable from setup_build.sh (x86_64-conda-linux-gnu-gcc).
# Always run `source /sandbox/setup_build.sh` before building.
# Make automatically picks up $CC from the environment.

CFLAGS = -std=c11 -Wall -Wextra -O2 -I include -I /sandbox/miniforge3/envs/build/include
LDFLAGS = -lm -lfftw3 -L/sandbox/miniforge3/envs/build/lib

SRCDIR  = src
INCDIR  = include
TESTDIR = tests
BUILDDIR = build

# Test programs (one per module)
TEST_NAMES = domain symmetry_ops basis initializers field engine analytic e2e python_compare

# Benchmark
BENCHMARK = benchmark

.PHONY: all clean test benchmark

all: $(TEST_NAMES)

# Create build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile test programs (they link directly against their .c files)
domain: $(TESTDIR)/test_domain.c $(SRCDIR)/domain.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

symmetry_ops: $(TESTDIR)/test_symmetry_ops.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

basis: $(TESTDIR)/test_basis.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c $(SRCDIR)/basis.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

initializers: $(TESTDIR)/test_initializers.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c $(SRCDIR)/basis.c $(SRCDIR)/initializers.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

field: $(TESTDIR)/test_field.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c $(SRCDIR)/basis.c $(SRCDIR)/initializers.c $(SRCDIR)/field.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

engine: $(TESTDIR)/test_engine.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c $(SRCDIR)/basis.c $(SRCDIR)/initializers.c $(SRCDIR)/field.c $(SRCDIR)/engine.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

analytic: $(TESTDIR)/test_analytic.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c $(SRCDIR)/basis.c $(SRCDIR)/initializers.c $(SRCDIR)/field.c $(SRCDIR)/engine.c $(SRCDIR)/analytic.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

e2e: $(TESTDIR)/test_e2e.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c $(SRCDIR)/basis.c $(SRCDIR)/initializers.c $(SRCDIR)/field.c $(SRCDIR)/engine.c $(SRCDIR)/analytic.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

python_compare: $(TESTDIR)/test_python_compare.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c $(SRCDIR)/basis.c $(SRCDIR)/initializers.c $(SRCDIR)/field.c $(SRCDIR)/engine.c $(SRCDIR)/analytic.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

# Benchmark
benchmark: $(TESTDIR)/benchmark.c $(SRCDIR)/domain.c $(SRCDIR)/symmetry_ops.c $(SRCDIR)/basis.c $(SRCDIR)/initializers.c $(SRCDIR)/field.c $(SRCDIR)/engine.c $(SRCDIR)/analytic.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

# Run all tests
test: all
	@for t in $(TEST_NAMES); do \
		echo "========== Running $$t =========="; \
		./$$t; \
		echo ""; \
	done
	@echo "========================================"
	@echo "All tests passed."

clean:
	rm -rf $(BUILDDIR) $(TEST_NAMES)
