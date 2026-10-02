#
# test/unit/bio/util/test_debruijn_graph.rb - Unit test for Bio::DeBruijnGraph
#
# Author::    BioRuby project
# License::   The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/util/debruijn_graph'
require 'bio/sequence'

module Bio
  class TestDeBruijnGraph < Test::Unit::TestCase
    def setup
      @klass = Bio::DeBruijnGraph
    end

    def build(k, reads)
      graph = @klass.new(k)
      reads.each { |read| graph.add(read) }
      graph
    end

    def test_new
      assert_equal(31, @klass.new(31).k)
    end

    def test_new_invalid_k
      assert_raise(ArgumentError) { @klass.new(1) }
      assert_raise(ArgumentError) { @klass.new(0) }
      assert_raise(ArgumentError) { @klass.new(-3) }
      assert_raise(ArgumentError) { @klass.new('31') }
    end

    def test_add_returns_self_for_chaining
      graph = @klass.new(3)
      graph << 'AGTCC' << 'AGTCC'
      assert_equal(2, graph.kmer_count('AGT'))
    end

    def test_add_accepts_sequence_objects_and_lowercase
      graph = @klass.new(3)
      graph.add(Bio::Sequence::NA.new('agtcc'))
      assert_equal(1, graph.kmer_count('AGT'))
    end

    def test_single_read_single_unitig
      unitigs = build(3, ['AGTCC']).unitigs
      assert_equal(1, unitigs.size)
      unitig = unitigs.first
      assert_equal('AGTCC', unitig.sequence)
      assert_equal([1, 1, 1], unitig.counts)
      assert_equal(1.0, unitig.coverage)
      assert_equal(3, unitig.kmer_count)
      assert_equal(5, unitig.length)
      assert_equal(false, unitig.circular?)
    end

    def test_coverage_counts_accumulate
      unitigs = build(3, ['AGTCC', 'GTCCA']).unitigs
      assert_equal(%w[AGTCCA], unitigs.collect(&:sequence))
      assert_equal([1, 2, 2, 1], unitigs[0].counts)
      assert_equal(6.0 / 4, unitigs[0].coverage)
      assert_equal(6, unitigs[0].total_count)
    end

    def test_kmer_count_is_canonical
      graph = build(3, ['AGTCC'])
      assert_equal(1, graph.kmer_count('AGT'))
      assert_equal(1, graph.kmer_count('ACT')) # reverse complement of AGT
      assert_equal(1, graph.kmer_count('TCC'))
      assert_equal(1, graph.kmer_count('GGA')) # reverse complement of TCC
    end

    def test_kmer_count_invalid_queries
      graph = build(3, ['AGTCC'])
      assert_equal(0, graph.kmer_count('AG'))
      assert_equal(0, graph.kmer_count('AGTC'))
      assert_equal(0, graph.kmer_count('ANG'))
      assert_equal(0, graph.kmer_count('AAA'))
    end

    def test_ambiguous_bases_break_windows
      graph = build(3, ['ACGANNNAGCA'])
      # no k-mer may span the N stretch
      assert_equal(4, graph.edge_count)
      assert_equal(%w[ACG AGC CGA GCA], graph.unitigs.collect(&:sequence))
    end

    def test_read_shorter_than_k_and_empty_reads
      graph = build(5, ['ACG', '', 'NNNNN'])
      assert_equal(0, graph.edge_count)
      assert_equal(0, graph.node_count)
      assert_equal([], graph.unitigs)
    end

    def test_reverse_complement_read_makes_no_new_contig
      forward_only = build(3, ['AGTCC']).unitigs
      both_strands = build(3, ['AGTCC', 'GGACT']).unitigs
      assert_equal(forward_only.collect(&:sequence),
                   both_strands.collect(&:sequence))
      assert_equal(['AGTCC'], both_strands.collect(&:sequence))
      assert_equal([2, 2, 2], both_strands.first.counts)
    end

    def test_junction_splits_unitigs
      unitigs = build(3, ['CAGT', 'CAGC']).unitigs
      assert_equal(%w[ACT AGC CAG], unitigs.collect(&:sequence))
      cag = unitigs.find { |u| u.sequence == 'CAG' }
      assert_equal([2], cag.counts)
    end

    def test_endpoint_orientations
      unitig = build(3, ['TTGA']).unitigs.first
      # reported in canonical orientation: min(TTGA, TCAA)
      assert_equal('TCAA', unitig.sequence)
      assert_equal('GA', unitig.left_node)
      assert_equal(:rev, unitig.left_orientation)
      assert_equal('AA', unitig.right_node)
      assert_equal(:fwd, unitig.right_orientation)
    end

    def test_isolated_cycle_has_stable_representation
      reads = %w[ACGT GTAC CGTA TACG]
      unitigs = build(3, reads).unitigs
      assert_equal(1, unitigs.size)
      unitig = unitigs.first
      assert_equal(true, unitig.circular?)
      assert_equal('ACGTAC', unitig.sequence)
      assert_equal([4, 4, 4, 4], unitig.counts)
      assert_equal(4.0, unitig.coverage)
      # rotated and shuffled reads give the same representation
      rotated = build(3, %w[TACG CGTA]).unitigs
      assert_equal(['ACGTAC'], rotated.collect(&:sequence))
      assert_equal([2, 2, 2, 2], rotated.first.counts)
    end

    def test_palindromic_kmer_orientation_is_stable
      # ACGT is a palindromic k-mer (k = 4, even)
      forward = build(4, ['AACGTA']).unitigs
      assert_equal(%w[AACG ACGT CGTA], forward.collect(&:sequence))
      # the reverse complement read must yield the identical contig set
      reverse = build(4, ['TACGTT']).unitigs
      assert_equal(forward.collect(&:sequence), reverse.collect(&:sequence))
      both = build(4, %w[AACGTA TACGTT])
      assert_equal(2, both.kmer_count('ACGT'))
      assert_equal(%w[AACG ACGT CGTA], both.unitigs.collect(&:sequence))
    end

    def test_shuffle_invariance
      reference = 'AGTCCAGTCTAGGCAACTGACCAGT'
      reads = []
      0.step(reference.length - 7, 3) do |i|
        reads << reference[i, 7]
      end
      reads += reads.collect { |r| Bio::DeBruijnGraph.reverse_complement(r) }
      graphs = [7, 42, 123].collect do |seed|
        build(4, reads.shuffle(random: Random.new(seed)))
      end
      expected_sequences = graphs.first.unitigs.collect(&:sequence)
      expected_counts = graphs.first.unitigs.collect(&:counts)
      graphs.each do |graph|
        assert_equal(expected_sequences, graph.unitigs.collect(&:sequence))
        assert_equal(expected_counts, graph.unitigs.collect(&:counts))
        assert_equal(graphs.first.kmer_counts, graph.kmer_counts)
      end
    end

    def test_node_and_edge_counts
      graph = build(3, ['AGTCC'])
      assert_equal(3, graph.edge_count)
      assert_equal(4, graph.node_count)
      assert_equal(%w[ACT GAC GGA].sort, graph.kmer_counts.keys.sort)
    end

    def test_reverse_complement_class_method
      assert_equal('GGACT', Bio::DeBruijnGraph.reverse_complement('AGTCC'))
      assert_equal('ACGT', Bio::DeBruijnGraph.reverse_complement('ACGT'))
    end

    def test_canonical_class_method
      assert_equal('ACT', Bio::DeBruijnGraph.canonical('ACT'))
      assert_equal('ACT', Bio::DeBruijnGraph.canonical('AGT'))
      assert_equal('ACGT', Bio::DeBruijnGraph.canonical('ACGT'))
    end
  end
end
