# frozen_string_literal: true

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/alignment'

module Bio
  module Alignment
    class TestNeiGojobori < Test::Unit::TestCase
      def test_calculates_nei_gojobori_distances_and_counts
        result = NeiGojobori.calculate('ATGCGACGA', 'ATGAGACGA')

        assert_equal(:ok, result.status)
        assert_equal([], result.status_reasons)
        assert_in_delta(6.6666666667, result.nonsynonymous_sites, 1e-9)
        assert_in_delta(2.3333333333, result.synonymous_sites, 1e-9)
        assert_in_delta(0.0, result.nonsynonymous_nucleotide_differences, 1e-12)
        assert_in_delta(1.0, result.synonymous_nucleotide_differences, 1e-12)
        assert_in_delta(0.0, result.dn, 1e-12)
        assert_in_delta(0.6354733953, result.ds, 1e-9)
        assert_in_delta(0.0, result.omega, 1e-12)
        assert_equal(3, result.comparable_codon_columns)
        assert_equal(0, result.excluded_columns)
      end

      def test_calculation_accepts_rna_and_a_codon_table_object
        table = Bio::CodonTable[1]
        result = NeiGojobori.new(codon_table: table).calculate('augcgacga',
                                                                'augagacga')

        assert_equal(:ok, result.status)
        assert_in_delta(0.6354733953, result.ds, 1e-9)
      end

      def test_multi_position_changes_are_averaged_over_valid_paths
        result = NeiGojobori.calculate('AAACGA', 'TTTCGA')

        assert_equal(:ok, result.status)
        assert_in_delta(4.3333333333, result.nonsynonymous_sites, 1e-9)
        assert_in_delta(1.6666666667, result.synonymous_sites, 1e-9)
        assert_in_delta(2.75, result.nonsynonymous_nucleotide_differences, 1e-12)
        assert_in_delta(0.25, result.synonymous_nucleotide_differences, 1e-12)
        assert_equal(2, result.comparable_codon_columns)
      end

      def test_invalid_replacement_paths_through_stops_are_ignored
        result = NeiGojobori.calculate('AAATTT', 'TACTTT')

        assert_equal(:zero_synonymous_distance, result.status)
        assert_in_delta(2.0, result.nonsynonymous_nucleotide_differences, 1e-12)
        assert_in_delta(0.0, result.synonymous_nucleotide_differences, 1e-12)
        assert_equal(0, result.excluded_columns)
      end

      def test_excludes_gap_ambiguous_and_stop_columns_by_default
        result = NeiGojobori.calculate('TTT---NNNTAGAAA',
                                       'TTTATGTAGTGGCGA')

        assert_equal(:ok, result.status)
        assert_equal(1, result.excluded_gap_columns)
        assert_equal(1, result.excluded_ambiguous_codon_columns)
        assert_equal(1, result.excluded_stop_codon_columns)
        assert_equal(3, result.excluded_columns)
        assert_equal(2, result.comparable_codon_columns)
      end

      def test_modes_can_raise_for_excluded_columns
        calculator = NeiGojobori.new(gap_mode: :raise,
                                    ambiguous_codon_mode: :raise,
                                    stop_codon_mode: :raise)

        assert_raise(ArgumentError) do
          calculator.calculate('---AAA', 'ATGAAA')
        end
        assert_raise(ArgumentError) do
          calculator.calculate('NNNAAA', 'ATGAAA')
        end
        assert_raise(ArgumentError) do
          calculator.calculate('TAAAAA', 'ATGAAA')
        end
      end

      def test_invalid_mode_is_rejected
        assert_raise(ArgumentError) do
          NeiGojobori.new(ambiguous_codon_mode: :include)
        end
      end

      def test_zero_synonymous_distance_has_an_explicit_status
        result = NeiGojobori.calculate('TGG', 'GGG')

        assert_equal(:zero_synonymous_distance, result.status)
        assert_equal(0.0, result.ds)
        assert_not_nil(result.dn)
        assert_nil(result.omega)
      end

      def test_no_synonymous_sites_has_an_explicit_status
        result = NeiGojobori.calculate('ATG', 'ATG')

        assert_equal(:undefined_synonymous_distance, result.status)
        assert_include(result.status_reasons, :no_synonymous_sites)
        assert_nil(result.ds)
        assert_nil(result.omega)
      end

      def test_jukes_cantor_out_of_domain_has_an_explicit_status
        result = NeiGojobori.calculate('TGG', 'GGC')

        assert_equal(:jukes_cantor_out_of_domain, result.status)
        assert_include(result.status_reasons, :jukes_cantor_out_of_domain)
        assert_not_nil(result.dn)
        assert_nil(result.omega)
      end

      def test_alignment_integration
        alignment = Alignment::SequenceArray['ATGCGACGA', 'ATGAGACGA']

        result = alignment.nei_gojobori_dnds

        assert_instance_of(NeiGojobori::Result, result)
        assert_equal(:ok, result.status)
        assert_in_delta(1.0, result.synonymous_nucleotide_differences, 1e-12)
      end

      def test_alignment_integration_uses_configured_gap_expression
        alignment = Alignment::SequenceArray['TTTNNN---', 'TTTATG###']
        alignment.gap_regexp = /#/

        result = alignment.nei_gojobori_dnds

        assert_equal(1, result.excluded_gap_columns)
        assert_equal(1, result.excluded_ambiguous_codon_columns)
      end

      def test_alignment_integration_accepts_direct_sequences
        alignment = Alignment::SequenceArray['ATGCGACGA', 'ATGAGACGA']

        result = alignment.nei_gojobori_dnds('ATGCGACGA', 'ATGAGACGA')

        assert_equal(:ok, result.status)
      end

      def test_original_alignment_integration
        alignment = Alignment::OriginalAlignment.new
        alignment.store('human', 'ATGCGACGA')
        alignment.store('mouse', 'ATGAGACGA')

        result = alignment.nei_gojobori_dnds('human', 'mouse')

        assert_equal(:ok, result.status)
        assert_in_delta(0.6354733953, result.ds, 1e-9)
      end

      def test_requires_equal_lengths_divisible_by_three
        assert_raise(ArgumentError) do
          NeiGojobori.calculate('ATGAA', 'ATGAAA')
        end
        assert_raise(ArgumentError) do
          NeiGojobori.calculate('ATGAAA', 'ATGAA')
        end
      end
    end
  end
end
