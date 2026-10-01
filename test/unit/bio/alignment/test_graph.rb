# frozen_string_literal: true

#
# test/unit/bio/alignment/test_graph.rb - unit tests for
# Bio::Alignment::Graph and Bio::Alignment::GraphAligner
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
require 'bio/alignment/graph'

module Bio
  # Independent pairwise aligner (Gotoh algorithm with affine gaps)
  # used as a reference to cross-check graph alignment scores.
  class TestGraphReferenceAligner
    NEGATIVE_INFINITY = -Float::INFINITY

    def initialize(match:, mismatch:, gap_open:, gap_extend:, mode:)
      @match = match
      @mismatch = mismatch
      @gap_open = gap_open
      @gap_extend = gap_extend
      @mode = mode
    end

    def score(str1, str2)
      @len1 = str1.length
      @len2 = str2.length
      @tables = Array.new(3) { Array.new(@len1 + 1) { Array.new(@len2 + 1, NEGATIVE_INFINITY) } }
      seed
      (0..@len1).each { |i| (0..@len2).each { |j| fill_cell(str1, str2, i, j) } }
      final_score
    end

    private

    def seed
      @tables[0][0][0] = 0
      return unless @mode == :global

      (1..@len1).each { |i| @tables[1][i][0] = @gap_open + (i - 1) * @gap_extend }
      (1..@len2).each { |j| @tables[2][0][j] = @gap_open + (j - 1) * @gap_extend }
    end

    def fill_cell(str1, str2, i, j)
      fill_insertion(i, j) if i.positive?
      fill_deletion(i, j) if j.positive?
      fill_match(str1, str2, i, j) if i.positive? && j.positive?
    end

    def fill_insertion(i, j)
      candidates = [@tables[0][i - 1][j] + @gap_open, @tables[1][i - 1][j] + @gap_extend]
      @tables[1][i][j] = floor(candidates.max)
    end

    def fill_deletion(i, j)
      candidates = [@tables[0][i][j - 1] + @gap_open, @tables[2][i][j - 1] + @gap_extend]
      @tables[2][i][j] = floor(candidates.max)
    end

    def fill_match(str1, str2, i, j)
      best = @tables.map { |table| table[i - 1][j - 1] }.max
      sub = str1[i - 1] == str2[j - 1] ? @match : @mismatch
      @tables[0][i][j] = floor(best + sub)
    end

    def floor(value)
      @mode == :local ? [value, 0].max : value
    end

    def final_score
      return @tables.map { |table| table[@len1][@len2] }.max if @mode == :global

      @tables.flat_map { |table| table.flatten }.max
    end
  end

  class TestAlignmentGraph < Test::Unit::TestCase
    def setup
      @options = { match: 2, mismatch: -1, gap_open: -3, gap_extend: -1 }
    end

    def chain_graph(string)
      graph = Bio::Alignment::Graph.new
      graph.add_path(string)
      graph
    end

    # A -> B -> D and A -> C -> D, where D is shared by both branches.
    def diamond_graph(char_b = 'C', char_c = 'G')
      graph = Bio::Alignment::Graph.new
      a = graph.add_node('A')
      b = graph.add_node(char_b)
      c = graph.add_node(char_c)
      d = graph.add_node('T')
      graph.add_edge(a, b)
      graph.add_edge(a, c)
      graph.add_edge(b, d)
      graph.add_edge(c, d)
      [graph, a, b, c, d]
    end

    def reference_score(str1, str2, mode)
      TestGraphReferenceAligner.new(mode: mode, **@options).score(str1, str2)
    end

    def assert_valid_node_chain(nodes)
      nodes.each_cons(2) do |pred, succ|
        assert_include(pred.successors, succ,
                       "traceback jumps along a nonexistent edge #{pred} -> #{succ}")
      end
    end

    def test_add_node_requires_single_character
      graph = Bio::Alignment::Graph.new
      assert_raise(ArgumentError) { graph.add_node('AC') }
      assert_raise(ArgumentError) { graph.add_node('') }
    end

    def test_add_edge_rejects_foreign_nodes
      graph = Bio::Alignment::Graph.new
      other = Bio::Alignment::Graph.new
      node = other.add_node('A')
      assert_raise(ArgumentError) { graph.add_edge(node, node) }
    end

    def test_each_path_expands_all_paths
      graph, = diamond_graph
      assert_equal(%w[ACT AGT], graph.path_strings.sort)
    end

    def test_cycle_is_rejected
      graph = Bio::Alignment::Graph.new
      x = graph.add_node('A')
      y = graph.add_node('C')
      graph.add_edge(x, y)
      graph.add_edge(y, x)
      assert_raise(Bio::Alignment::Graph::CyclicError) { graph.topological_order }
      assert_raise(Bio::Alignment::Graph::CyclicError) { graph.each_path.to_a }
      assert_raise(Bio::Alignment::Graph::CyclicError) do
        Bio::Alignment::GraphAligner.new.align('A', graph)
      end
    end

    def test_invalid_mode_is_rejected
      assert_raise(ArgumentError) { Bio::Alignment::GraphAligner.new(mode: :semiglobal) }
    end

    def test_empty_graph_is_rejected
      aligner = Bio::Alignment::GraphAligner.new
      assert_raise(ArgumentError) { aligner.align('ACGT', Bio::Alignment::Graph.new) }
    end

    def test_global_chain_matches_reference
      graph = chain_graph('ACGTAC')
      result = Bio::Alignment::GraphAligner.new(mode: :global, **@options).align('AGTAC', graph)
      assert_equal(reference_score('ACGTAC', 'AGTAC', :global), result.score)
      assert_equal('ACGTAC', result.aligned_graph)
      assert_equal('A-GTAC', result.aligned_sequence)
    end

    def test_affine_gap_open_charged_once_across_nodes
      graph = chain_graph('AAAAAA')
      result = Bio::Alignment::GraphAligner.new(mode: :global, **@options).align('AA', graph)
      # 2 matches (2 * 2) and one continuous 4-node deletion (-3 + 3 * -1)
      assert_equal(-2, result.score)
      assert_equal(reference_score('AAAAAA', 'AA', :global), result.score)
      deleted = result.pairs.select { |_, index| index.nil? }
      assert_equal(4, deleted.size)
    end

    def test_insertion_against_graph
      graph = chain_graph('AC')
      result = Bio::Alignment::GraphAligner.new(mode: :global, **@options).align('AAC', graph)
      assert_equal(reference_score('AC', 'AAC', :global), result.score)
      assert_equal(1, result.pairs.count { |node, _| node.nil? })
    end

    def test_local_mode
      graph = chain_graph('ACGT')
      result = Bio::Alignment::GraphAligner.new(mode: :local, **@options).align('TTACGTAA', graph)
      assert_equal(8, result.score)
      assert_equal('ACGT', result.aligned_sequence)
      assert_equal(2, result.sequence_begin)
      assert_equal(6, result.sequence_end)
    end

    def test_local_mode_no_match_gives_empty_alignment
      graph = chain_graph('AAAA')
      result = Bio::Alignment::GraphAligner.new(mode: :local, **@options).align('TTTT', graph)
      assert_equal(0, result.score)
      assert_equal([], result.pairs)
    end

    def test_shared_node_is_not_spliced_into_other_branch
      # Both branches carry 'G', but there is no edge B -> C: the
      # traceback must never join the two branches by matching
      # characters.
      graph, a, b, = diamond_graph('G', 'G')
      result = Bio::Alignment::GraphAligner.new(mode: :global, **@options).align('AGGT', graph)
      # A spliced path A -> B -> C -> D ("AGGT") would score 8; the best
      # real path scores 2 * 3 + gap_open = 3.
      assert_equal(3, result.score)
      assert_equal([a, b, graph.nodes[3]], result.nodes)
      assert_valid_node_chain(result.nodes)
    end

    def test_tie_break_is_deterministic
      graph, a, b, = diamond_graph('G', 'G')
      aligner = Bio::Alignment::GraphAligner.new(mode: :global, **@options)
      first = aligner.align('AGT', graph)
      second = aligner.align('AGT', graph)
      assert_equal(first.pairs.map { |n, i| [n&.index, i] },
                   second.pairs.map { |n, i| [n&.index, i] })
      # the branch inserted first (b) wins the tie
      assert_equal([a, b, graph.nodes[3]], first.nodes)
    end

    def test_substitution_hash
      graph = chain_graph('AG')
      substitution = { %w[A G] => 3 }
      aligner = Bio::Alignment::GraphAligner.new(mode: :global, **@options,
                                                 substitution: substitution)
      assert_equal(2 + 3, aligner.align('AA', graph).score)
    end

    def test_random_dags_match_best_path_pairwise_score
      rng = Random.new(20261001)
      5.times do
        graph = random_dag(rng)
        3.times do
          query = random_sequence(rng, 6)
          %i[global local].each do |mode|
            assert_graph_score_equals_best_path(graph, query, mode)
          end
        end
      end
    end

    def random_dag(rng)
      graph = Bio::Alignment::Graph.new
      nodes = Array.new(7) { graph.add_node(%w[A C][rng.rand(2)]) }
      nodes.combination(2) do |u, v|
        graph.add_edge(u, v) if rng.rand < 0.35
      end
      graph
    end

    def random_sequence(rng, length)
      Array.new(length) { %w[A C][rng.rand(2)] }.join
    end

    def assert_graph_score_equals_best_path(graph, query, mode)
      aligner = Bio::Alignment::GraphAligner.new(mode: mode, **@options)
      result = aligner.align(query, graph)
      expected = graph.path_strings.map { |path| reference_score(path, query, mode) }.max
      assert_equal(expected, result.score,
                   "mode=#{mode} query=#{query} paths=#{graph.path_strings.inspect}")
      assert_valid_node_chain(result.nodes)
    end
  end
end
