# frozen_string_literal: true
#
# test/unit/bio/util/test_nei_gojobori.rb - Unit test for Bio::NeiGojobori
#
# Copyright::  Copyright (C) 2026 BioRuby project
# License::    The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/alignment'
require 'bio/util/nei_gojobori'

module Bio
  class TestNeiGojobori < Test::Unit::TestCase
    def setup
      @ng = Bio::NeiGojobori.new
    end

    def test_identical_sequences
      result = @ng.run('gctgctgct', 'gctgctgct')
      assert_equal(3, result.compared_columns)
      assert_equal(0, result.excluded_columns)
      assert_in_delta(3.0, result.synonymous_sites)
      assert_in_delta(6.0, result.nonsynonymous_sites)
      assert_in_delta(0.0, result.synonymous_differences)
      assert_in_delta(0.0, result.nonsynonymous_differences)
      assert_in_delta(0.0, result.dN)
      assert_in_delta(0.0, result.dS)
      assert_equal(:zero_ds, result.status)
      assert_nil(result.omega)
      assert_false(result.ok?)
    end

    def test_single_synonymous_difference
      # 10 identical gct columns + ttt/ttc (Phe/Phe, 1 synonymous change)
      seq1 = 'gct' * 10 + 'ttt'
      seq2 = 'gct' * 10 + 'ttc'
      result = @ng.run(seq1, seq2)
      assert_equal(11, result.compared_columns)
      assert_in_delta(31.0 / 3, result.synonymous_sites)
      assert_in_delta(68.0 / 3, result.nonsynonymous_sites)
      assert_in_delta(1.0, result.synonymous_differences)
      assert_in_delta(0.0, result.nonsynonymous_differences)
      assert_in_delta(3.0 / 31, result.pS)
      assert_in_delta(0.0, result.pN)
      expected_ds = -0.75 * Math.log(27.0 / 31)
      assert_in_delta(expected_ds, result.dS)
      assert_in_delta(0.0, result.dN)
      assert_equal(:ok, result.status)
      assert(result.ok?)
      assert_in_delta(0.0, result.omega)
      assert_in_delta(0.0, result.dN_dS)
    end

    def test_single_nonsynonymous_difference
      # 10 identical gct columns + gct/ggt (Ala/Gly, 1 nonsynonymous change)
      seq1 = 'gct' * 10 + 'gct'
      seq2 = 'gct' * 10 + 'ggt'
      result = @ng.run(seq1, seq2)
      assert_in_delta(1.0, result.nonsynonymous_differences)
      assert_in_delta(0.0, result.synonymous_differences)
      expected_dn = -0.75 * Math.log(31.0 / 33)
      assert_in_delta(expected_dn, result.dN)
      assert_in_delta(0.0, result.dS)
      assert_equal(:zero_ds, result.status)
      assert_nil(result.omega)
    end

    def test_omega
      # 10 identical gct + ttt/ttc (1 syn) + gct/ggt (1 nonsyn)
      seq1 = 'gct' * 10 + 'ttt' + 'gct'
      seq2 = 'gct' * 10 + 'ttc' + 'ggt'
      result = @ng.run(seq1, seq2)
      assert_in_delta(34.0 / 3, result.synonymous_sites)
      assert_in_delta(74.0 / 3, result.nonsynonymous_sites)
      assert_in_delta(3.0 / 34, result.pS)
      assert_in_delta(3.0 / 74, result.pN)
      expected_ds = -0.75 * Math.log(15.0 / 17)
      expected_dn = -0.75 * Math.log(35.0 / 37)
      assert_in_delta(expected_ds, result.dS)
      assert_in_delta(expected_dn, result.dN)
      assert_equal(:ok, result.status)
      assert_in_delta(expected_dn / expected_ds, result.omega)
    end

    def test_two_differences_averaged_over_paths
      # tta (Leu) vs ctt (Leu): one path is fully synonymous
      # (tta -> cta -> ctt), the other fully nonsynonymous
      # (tta -> ttt -> ctt); the average must be 1/1, not 2/0 or 0/2.
      result = @ng.run('tta', 'ctt')
      assert_in_delta(1.0, result.synonymous_differences)
      assert_in_delta(1.0, result.nonsynonymous_differences)
      assert_in_delta(5.0 / 6, result.synonymous_sites)
      assert_in_delta(13.0 / 6, result.nonsynonymous_sites)
      # pS = 1 / (5/6) = 1.2 is outside the Jukes-Cantor domain
      assert_equal(:saturated, result.status)
      assert(result.saturated?)
      assert_nil(result.dS)
      assert_nil(result.omega)
    end

    def test_three_differences_skip_paths_through_stop
      # tgg (Trp) vs cta (Leu): 2 of the 6 mutational orders pass
      # through the tga stop codon and must be discarded; averaging
      # over the 4 valid paths gives 3/2 synonymous and 3/2
      # nonsynonymous differences.
      result = @ng.run('tgg', 'cta')
      assert_in_delta(1.5, result.synonymous_differences)
      assert_in_delta(1.5, result.nonsynonymous_differences)
      assert_equal(1, result.compared_columns)
      assert_equal(0, result.excluded_columns)
    end

    def test_gap_column_excluded_by_default
      result = @ng.run('gct---', 'gctgct')
      assert_equal(1, result.compared_columns)
      assert_equal(1, result.excluded_columns)
    end

    def test_gap_column_raise
      ng = Bio::NeiGojobori.new(gap: :raise)
      assert_raise(Bio::NeiGojobori::Error) do
        ng.run('gct---', 'gctgct')
      end
    end

    def test_ambiguous_column_excluded_by_default
      result = @ng.run('gctgct', 'gctgcn')
      assert_equal(1, result.compared_columns)
      assert_equal(1, result.excluded_columns)
    end

    def test_ambiguous_column_raise
      ng = Bio::NeiGojobori.new(ambiguous: :raise)
      assert_raise(Bio::NeiGojobori::Error) do
        ng.run('gctgct', 'gctgcn')
      end
    end

    def test_stop_column_excluded_by_default
      result = @ng.run('gctgct', 'gcttaa')
      assert_equal(1, result.compared_columns)
      assert_equal(1, result.excluded_columns)
    end

    def test_stop_column_raise
      ng = Bio::NeiGojobori.new(stop: :raise)
      assert_raise(Bio::NeiGojobori::Error) do
        ng.run('gctgct', 'gcttaa')
      end
    end

    def test_all_columns_excluded_gives_no_sites
      result = @ng.run('---', 'gct')
      assert_equal(0, result.compared_columns)
      assert_equal(1, result.excluded_columns)
      assert_equal(:no_sites, result.status)
      assert_nil(result.dN)
      assert_nil(result.dS)
      assert_nil(result.omega)
    end

    def test_custom_codon_table
      # tga is a stop in the standard code but Trp in the vertebrate
      # mitochondrial code (table 2).
      standard = @ng.run('tgg', 'tga')
      assert_equal(0, standard.compared_columns)
      assert_equal(1, standard.excluded_columns)

      mito = Bio::NeiGojobori.new(codon_table: Bio::CodonTable[2])
      result = mito.run('tgg', 'tga')
      assert_equal(1, result.compared_columns)
      assert_in_delta(1.0, result.synonymous_differences)
      assert_in_delta(0.0, result.nonsynonymous_differences)
    end

    def test_run_alignment
      alignment = Bio::Alignment.new(['gctgct', 'gctggt'])
      result = @ng.run_alignment(alignment)
      assert_equal(2, result.compared_columns)
      assert_in_delta(1.0, result.nonsynonymous_differences)
    end

    def test_run_alignment_class_method
      alignment = Bio::Alignment.new(['gctgct', 'gctggt'])
      result = Bio::NeiGojobori.run_alignment(alignment)
      assert_equal(2, result.compared_columns)
    end

    def test_run_alignment_requires_two_sequences
      alignment = Bio::Alignment.new(%w[gct gct gct])
      assert_raise(Bio::NeiGojobori::Error) do
        @ng.run_alignment(alignment)
      end
    end

    def test_class_method_run
      result = Bio::NeiGojobori.run('gct', 'ggt')
      assert_in_delta(1.0, result.nonsynonymous_differences)
    end

    def test_unequal_lengths_raise
      assert_raise(Bio::NeiGojobori::Error) do
        @ng.run('gctgct', 'gct')
      end
    end

    def test_length_not_multiple_of_three_raises
      assert_raise(Bio::NeiGojobori::Error) do
        @ng.run('gctg', 'gctg')
      end
    end

    def test_invalid_policy_raises_argument_error
      assert_raise(ArgumentError) do
        Bio::NeiGojobori.new(gap: :ignore)
      end
    end

    def test_accepts_bio_sequence_na_and_uppercase
      seq1 = Bio::Sequence::NA.new('GCTGCT')
      seq2 = Bio::Sequence::NA.new('GCTGGT')
      result = @ng.run(seq1, seq2)
      assert_in_delta(1.0, result.nonsynonymous_differences)
    end
  end
end
