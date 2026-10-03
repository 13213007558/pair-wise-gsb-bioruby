# coding: US-ASCII
# frozen_string_literal: true

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio'

module Bio
  class TestVCFReader < Test::Unit::TestCase
    def setup
      @data = <<~END_OF_VCF
        ##fileformat=VCFv4.2
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo
        chr1\t10\t.\tA\tT,C\t30\tPASS\tAF=0.25,0.5\tGT:PL\t0|1:0,10,20,10,20,30\t2/1:30,20,10,20,0,10
        chr1\t12\tdb\tG\t<DEL>\t.\t.\t.\tGT\t.\t0/1
      END_OF_VCF
      @vcf = Bio::VCF.new(@data)
    end

    def test_parses_header_and_records
      assert_equal(2, @vcf.records.length)
      assert_equal(%w[one two], @vcf.header.sample_names)
      assert_equal('A', @vcf.header.info['AF']['Number'])
    end

    def test_parses_info_format_and_samples
      record = @vcf.records.first
      assert_equal(%w[T C], record.alt)
      assert_equal('0.25,0.5', record.info['AF'])
      assert_equal(%w[GT PL], record.format)
      assert_equal('0|1', record.samples['one']['GT'])
      assert_equal('30,20,10,20,0,10', record.samples['two']['PL'])
    end

    def test_round_trips_vcf
      assert_equal(@data, @vcf.to_s)
    end

    def test_detects_symbolic_record
      assert_equal(false, @vcf.records.first.symbolic?)
      assert_equal(true, @vcf.records.last.symbolic?)
    end

    def test_records_without_header_round_trip_without_synthetic_header
      line = "chr1\t10\t.\tA\tT\t.\tPASS\t.\n"
      assert_equal(line, Bio::VCF.new(line).to_s)
    end
  end

  class TestVCFNormalize < Test::Unit::TestCase
    def test_left_aligns_deletion_in_repeat
      record = record("chr1\t2\t.\tTT\tT", '.', 'GT', '1/1')
      normalized = Bio::VCF::Normalizer.new({ 'chr1' => 'TTT' })
                                           .normalize(record)

      assert_equal(1, normalized.pos)
      assert_equal('TT', normalized.ref)
      assert_equal(['T'], normalized.alt)
      assert_equal(:normalized, normalized.normalization_status)
    end

    def test_equivalent_deletion_representations_normalize_equally
      first = record("chr1\t1\t.\tAT\tA", '.', 'GT', '0/1')
      second = record("chr1\t2\t.\tTT\tT", '.', 'GT', '0/1')
      reference = { 'chr1' => 'ATTT' }

      assert_equal(normalized(first, reference),
                   normalized(second, reference))
    end

    def test_left_aligns_insertion_in_repeat
      record = record("chr1\t3\t.\tT\tTT", '.', 'GT', '0/1')
      normalized = Bio::VCF::Normalizer.new({ 'chr1' => 'TTT' })
                                           .normalize(record)

      assert_equal(1, normalized.pos)
      assert_equal('T', normalized.ref)
      assert_equal(['TT'], normalized.alt)
    end

    def test_minimizes_substitution
      record = record("chr1\t10\t.\tGAT\tGCT", '.', 'GT', '0/1')
      normalized = Bio::VCF::Normalizer.new({ 'chr1' => 'NGAT' },
                                           sequence_start: 9).normalize(record)

      assert_equal(11, normalized.pos)
      assert_equal('A', normalized.ref)
      assert_equal(['C'], normalized.alt)
    end

    def test_leaves_symbolic_allele_unprocessed
      symbolic = record("chr1\t10\t.\tG\t<DEL>", '.', 'GT', '0/1')
      normalized = Bio::VCF::Normalizer.new({ 'chr1' => 'GAT' })
                                           .normalize(symbolic)

      assert_equal(10, normalized.pos)
      assert_equal('G', normalized.ref)
      assert_equal(['<DEL>'], normalized.alt)
      assert_equal(:unprocessed, normalized.normalization_status)
    end

    def test_split_reorders_info_and_sample_allele_fields
      data = <<~END_OF_VCF
        ##INFO=<ID=AA,Number=A,Type=Integer,Description="A">
        ##INFO=<ID=RR,Number=R,Type=Integer,Description="R">
        ##FORMAT=<ID=XX,Number=G,Type=Integer,Description="G">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
        chr1\t2\t.\tT\tA,C\t.\tPASS\tAA=10,20;RR=3,4,5;FLAG\tGT:AD:XX:PL\t0|2:3,4,5:0,1,2,3,4,5:9,8,7,6,5,0
      END_OF_VCF
      parsed = Bio::VCF.new(data)
      reference = { 'chr1' => 'TTAC' }
      records = Bio::VCF::Normalizer.new(reference).normalize(parsed.records.first,
                                                            split: true)

      assert_equal(2, records.length)
      assert_equal(['A'], records.first.alt)
      assert_equal(['C'], records.last.alt)
      assert_equal('10', records.first.info['AA'])
      assert_equal('20', records.last.info['AA'])
      assert_equal('3,4', records.first.info['RR'])
      assert_equal('3,5', records.last.info['RR'])
      assert_equal(true, records.first.info['FLAG'])

      first_sample = records.first.samples['sample']
      last_sample = records.last.samples['sample']
      assert_equal('0|.', first_sample['GT'])
      assert_equal('0|1', last_sample['GT'])
      assert_equal('3,4', first_sample['AD'])
      assert_equal('3,5', last_sample['AD'])
      assert_equal('0,1,3', first_sample['XX'])
      assert_equal('0,2,5', last_sample['XX'])
      assert_equal('3,2,0', first_sample['PL'])
      assert_equal('9,7,0', last_sample['PL'])
    end

    private

    def record(pos_id_ref_alt, info, format_key, genotype)
      chrom, pos, id, ref, alt = pos_id_ref_alt.split("\t")
      line = [chrom, pos, id, ref, alt, '.', '.', info, format_key,
              genotype].join("\t") + "\n"
      header = Bio::VCF::Header.new
      header.sample_names = ['sample']
      Bio::VCF::Record.new(line, header)
    end

    def normalized(record, reference)
      Bio::VCF::Normalizer.new(reference).normalize(record).to_s
    end
  end
end
