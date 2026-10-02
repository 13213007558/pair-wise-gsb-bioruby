require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 2,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio'
require 'bio/assembly'

module Bio
  class TestAssemblyUnitig < Test::Unit::TestCase
    def test_builds_maximal_linear_unitig_with_coverage
      graph = Assembly::DeBruijnGraph.new(3, %w[ATGCA ATGCA])
      unitig = graph.unitigs.first

      assert_equal(['atgca'], graph.unitigs.map(&:sequence))
      assert_equal(false, unitig.circular?)
      assert_equal([2, 4, 4], unitig.kmer_counts)
      assert_equal([2, 2], unitig.edge_counts)
      assert_equal('atg', unitig.start_kmer)
      assert_equal('gca', unitig.end_kmer)
      assert_equal('+', unitig.start_orientation)
      assert_equal('+', unitig.end_orientation)
      assert_equal(2.0, unitig.coverage)
    end

    def test_reverse_complement_reads_do_not_make_extra_contig
      graph = Assembly::DeBruijnGraph.new(3, %w[ATGCA TGCAT])
      unitig = graph.unitigs.first

      assert_equal(['atgca'], graph.unitigs.map(&:sequence))
      assert_equal([2, 4, 4], unitig.kmer_counts)
      assert_equal([2, 2], unitig.edge_counts)
    end

    def test_ambiguous_bases_cut_windows
      graph = Assembly::DeBruijnGraph.new(3, ['NATGCxTGCAN ATGCA'])
      unitig = graph.unitigs.first

      assert_equal(['atgca'], graph.unitigs.map(&:sequence))
      assert_equal([2, 5, 5], unitig.kmer_counts)
      assert_equal([2, 2], unitig.edge_counts)
    end

    def test_shorter_windows_are_ignored_and_short_reads_are_isolated
      graph = Assembly::DeBruijnGraph.new(3, ['ac', 'ATG'])
      unitig = graph.unitigs.first

      assert_equal(['atg'], graph.unitigs.map(&:sequence))
      assert_equal([1], unitig.kmer_counts)
      assert_equal([], unitig.edge_counts)
      assert_equal(1.0, unitig.coverage)
    end

    def test_branching_reads_produce_separate_unitigs
      graph = Assembly::DeBruijnGraph.new(3, ['ATGCA', 'ATGTA'])
      result = graph.unitigs.map { |unitig| [unitig.sequence, unitig.coverage] }

      assert_equal([['atgca', 1.0], ['atgta', 1.0]], result)
    end

    def test_input_order_does_not_change_sequence_set_or_coverage
      reads = ['ATGCA', 'CATGC', 'TGCAT', 'GCATG']
      expected = Assembly.unitigs(reads, 3)
      shuffled = Assembly.unitigs(reads.reverse, 3)

      assert_equal(expected.map(&:sequence), shuffled.map(&:sequence))
      assert_equal(expected.map(&:kmer_counts), shuffled.map(&:kmer_counts))
      assert_equal(expected.map(&:edge_counts), shuffled.map(&:edge_counts))
    end

    def test_isolated_cycle_has_stable_representation
      graph = Assembly::DeBruijnGraph.new(2, ['AAGTCCAA'])
      unitig = graph.unitigs.first

      assert_equal(1, graph.unitigs.size)
      assert_equal(true, unitig.circular?)
      assert_equal('aagtcca', unitig.sequence)
      assert_equal([2, 1, 1, 1, 1, 1], unitig.kmer_counts)
      assert_equal([1] * 6, unitig.edge_counts)
      assert_equal('aa', unitig.start_kmer)
      assert_equal('aa', unitig.end_kmer)
    end

    def test_one_base_isolated_cycle_does_not_duplicate_kmer_coverage
      graph = Assembly::DeBruijnGraph.new(1, ['AA'])
      unitig = graph.unitigs.first

      assert_equal(['a'], graph.unitigs.map(&:sequence))
      assert_equal(true, unitig.circular?)
      assert_equal([2], unitig.kmer_counts)
      assert_equal([1], unitig.edge_counts)
    end

    def test_even_k_palindrome_orientation_is_deterministic
      graph = Assembly::DeBruijnGraph.new(2, ['AT', 'at', 'At'])
      unitig = graph.unitigs.first

      assert_equal(['at'], graph.unitigs.map(&:sequence))
      assert_equal(false, unitig.circular?)
      assert_equal([3], unitig.kmer_counts)
      assert_equal([], unitig.edge_counts)
      assert_equal('+', unitig.start_orientation)
      assert_equal('+', unitig.end_orientation)
    end

    def test_even_k_palindromic_path_uses_stable_endpoint_orientation
      forward = Assembly::DeBruijnGraph.new(2, ['CATG']).unitigs
      reverse = Assembly::DeBruijnGraph.new(2, ['CATG']).unitigs

      assert_equal(1, forward.size)
      assert_equal(1, reverse.size)
      assert_equal('atg', forward.first.sequence)
      assert_equal('atg', reverse.first.sequence)
      assert_equal('at', forward.first.start_kmer)
      assert_equal('tg', forward.first.end_kmer)
      assert_equal(
        forward.map(&:unitig_key), reverse.map(&:unitig_key)
      )
    end

    def test_rna_is_normalized
      graph = Assembly::DeBruijnGraph.new(3, ['AUGCA'])

      assert_equal(['atgca'], graph.unitigs.map(&:sequence))
    end

    def test_invalid_k_is_rejected
      assert_raise(ArgumentError) { Assembly::DeBruijnGraph.new(0) }
    end

    def test_empty_input_returns_no_unitigs
      assert_equal([], Assembly.unitigs([], 3))
    end
  end
end
