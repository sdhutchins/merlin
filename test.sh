#!/bin/bash
# MERLIN Test Script
# Tests all examples mentioned in README.md
set -e  # Exit on any error

echo "MERLIN Test Script - Testing README Examples"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to test a command
test_command() {
    local description="$1"
    local command="$2"
    local expected_exit="$3"
    
    echo -e "\n${YELLOW}Testing: $description${NC}"
    echo "Command: $command"
    
    if eval "$command" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASSED${NC}"
        ((passed++))
        return 0
    else
        echo -e "${RED}✗ FAILED${NC}"
        ((failed++))
        return 1
    fi
}

# Function to test file creation
test_file_created() {
    local file="$1"
    local description="$2"
    
    if [ -f "$file" ]; then
        echo -e "${GREEN}✓ PASSED: $description - File created: $file${NC}"
        ((passed++))
        return 0
    else
        echo -e "${RED}✗ FAILED: $description - File not created: $file${NC}"
        ((failed++))
        return 1
    fi
}

# Check if executables exist
echo -e "\n${YELLOW}Checking executables...${NC}"
if [ ! -f "./executables/merlin" ]; then
    echo -e "${RED}✗ merlin executable not found. Run 'make' first.${NC}"
    exit 1
fi

if [ ! -f "./executables/pedstats" ]; then
    echo -e "${RED}✗ pedstats executable not found. Run 'make' first.${NC}"
    exit 1
fi

echo -e "${GREEN}✓ All executables found${NC}"

# Test counter
passed=0
failed=0

# Test 1: Basic pedstats example
echo -e "\n${YELLOW}=== Testing Pedstats Examples ===${NC}"

test_command "Basic pedigree statistics" \
    "./executables/pedstats -d examples/basic2.dat -p examples/basic2.ped"

# Test 2: Hardy-Weinberg testing
test_command "Hardy-Weinberg testing" \
    "./executables/pedstats -d examples/basic2.dat -p examples/basic2.ped --hardyWeinberg"

# Test 3: PDF generation
test_command "PDF generation" \
    "./executables/pedstats -d examples/basic2.dat -p examples/basic2.ped --pdf"

# Check if PDF was created
test_file_created "pedstats.pdf" "PDF output file"

# Test 4: Large dataset analysis
test_command "Large dataset analysis" \
    "./executables/pedstats -d examples/assoc.dat -p examples/assoc.ped --verbose"

# Test 5: Basic merlin analysis
echo -e "\n${YELLOW}=== Testing Merlin Examples ===${NC}"

test_command "Basic merlin error checking" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --error"

test_command "Basic merlin information content" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --info"

# Test 6: Cluster-based SNP analysis
test_command "Cluster-based SNP analysis" \
    "./executables/merlin -d examples/snp-scan.dat -p examples/snp-scan.ped -m examples/snp-scan.map --clusters examples/snp-scan.clusters"

# Test 7: Parametric linkage analysis
test_command "Parametric linkage analysis" \
    "./executables/merlin -d examples/parametric.dat -p examples/parametric.ped -m examples/parametric.map --model examples/parametric.model --freq examples/parametric.freq"

# Test 8: Additional merlin examples from README
echo -e "\n${YELLOW}=== Testing Additional Merlin Features ===${NC}"

test_command "IBD and kinship matrices" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --ibd --kinship"

test_command "Non-parametric linkage (pairs)" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --pairs"

test_command "Non-parametric linkage (npl)" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --npl"

test_command "QTL analysis" \
    "./executables/merlin -d examples/assoc.dat -p examples/assoc.ped -m examples/assoc.map --qtl"

test_command "Variance components analysis" \
    "./executables/merlin -d examples/assoc.dat -p examples/assoc.ped -m examples/assoc.map --vc"

test_command "Best inheritance vectors" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --best"

test_command "Sample inheritance vectors" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --sample"

# Test 9: Memory management examples
echo -e "\n${YELLOW}=== Testing Memory Management ===${NC}"

test_command "Memory limit (256MB)" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --megabytes:256"

test_command "Complexity limit (16 bits)" \
    "./executables/merlin -d examples/basic2.dat -p examples/basic2.ped -m examples/basic2.map --bits:16"

# Summary
echo -e "\n${YELLOW}=========================================="
echo "TEST SUMMARY"
echo "==========================================${NC}"
echo -e "${GREEN}Passed: $passed${NC}"
echo -e "${RED}Failed: $failed${NC}"
echo -e "Total: $((passed + failed))"

if [ $failed -eq 0 ]; then
    echo -e "\n${GREEN}🎉 All tests passed!${NC}"
    exit 0
else
    echo -e "\n${RED}❌ Some tests failed. Check the output above.${NC}"
    exit 1
fi 