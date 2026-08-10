#!/bin/bash
# Test suite to check RUFUS working with DAC 

# Fill in
ST002_1D_BCM_300x=/mnt/s3/test_files/ST002-1D_BCM_ill_DS300x0.bam
ST002_1D_BCM_300x_index=/mnt/s3/test_files/ST002-1D_BCM_ill_DS300x0.bam.bai
REF_DIR=/mnt/s3/references/
HASH_DIR=/mnt/s3/rufus_resources/
REGION_DIR=./region_files/
PIPELINE_TEST=./pipeline_test.sh

TEST_MOUNT="-v $ST002_1D_BCM_300x:$ST002_1D_BCM_300x \
-v $ST002_1D_BCM_300x_index:$ST002_1D_BCM_300x_index \
-v $REF_DIR:$REF_DIR \
-v $HASH_DIR:$HASH_DIR"

stop_container() {
  docker stop rufus-worker
}
trap 'stop_container' EXIT

export RUFUS_ROOT=/opt/RUFUS/

# Start container
USER_SPEC="$(id -u):$(id -g)"
CONTAINER_ID=$(docker run -d --rm --name rufus-worker \
  -u "${USER_SPEC}" \
  -v ".:/host" \
  --cap-add SYS_ADMIN \
  --device /dev/fuse \
  $input_mount_clause \
  $TEST_MOUNT \
  rufus:dac \
  tail -f /dev/null)

# Test 1
mkdir -p dac_single_region
bash "$PIPELINE_TEST" -s "$ST002_1D_BCM_300x" -r "$REF_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" -f "${REGION_DIR}/single_region.tsv" -i 1 -d "$HASH_DIR" -p "dac_single_region" -o "test_1"

# Test 2
mkdir -p dac_four_regions
bash "$PIPELINE_TEST" -s "$ST002_1D_BCM_300x" -r "$REF_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" -f "$REGION_DIR/four_regions.tsv" -i 1 -d "$HASH_DIR" -p "dac_four_regions" -o "test_2"

# Test 3
mkdir -p dac_two_shards_four_regions
bash "$PIPELINE_TEST" -s "$ST002_1D_BCM_300x" -r "$REF_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" -f "$REGION_DIR/two_shards_four_regions.tsv" -i 1 -d "$HASH_DIR" -p "dac_two_shards_four_regions" -o "test_3"
bash "$PIPELINE_TEST" -s "$ST002_1D_BCM_300x" -r "$REF_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" -f "$REGION_DIR/two_shards_four_regions.tsv" -i 2 -d "$HASH_DIR" -p "dac_two_shards_four_regions" -o "test_3"

# Test 4
mkdir -p dac_ten_shards_100_regions
for i in {1..10}; do
	bash "$PIPELINE_TEST" -s "$ST002_1D_BCM_300x" -r "$REF_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" -f "$REGION_DIR/ten_shards_100_regions.tsv" -i $i -d "$HASH_DIR" -p "dac_ten_shards_100_regions" -o "test_4"
done

# Count results
TEST_OUT="test_results.out"
echo -n "test_1.vcf.gz: "
bcftools view -H "dac_single_region/test_1.vcf.gz" | wc -l >> "$TEST_OUT"

echo -n "test_2.vcf.gz: "
bcftools view -H "dac_four_regions/test_2.vcf.gz" | wc -l >> "$TEST_OUT"

echo -n "test_3.vcf.gz: "
bcftools view -H "dac_two_shards_four_regions/test_3.vcf.gz" | wc -l >> "$TEST_OUT"

echo -n "test_4.vcf.gz: "
bcftools view -H "dac_ten_shards_100_regions/test_4.vcf.gz" | wc -l >> "$TEST_OUT"

# Check same data different shard numbers produce identical results
bcftools isec -p isecs -c none dac_four_regions/test_2.vcf.gz dac_two_shards_four_regions/test_3.vcf.gz
unique_single_shard_count=$(bcftools view -H 0000.vcf | wc -l)
unique_two_shard_count=$(bcftools view -H 0001.vcf | wc -l)
if [[ "$unique_single_shard_count" -eq 0 ]] && [[ "$unique_two_shard_count" -eq 0 ]]; then
	echo "Intersection test passed" >> "$TEST_OUT"
else
	echo "Intersection test failed: using two shards for the same regions does not equal using a single shard" >> "$TEST_OUT"
fi

echo "Tests completed, see $TEST_OUT for results."