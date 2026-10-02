#
# test/unit/bio/test_dcj.rb - Unit test for Bio::DCJ
#
# Copyright:: Copyright (C) 2026 BioRuby Project
# License:: The Ruby License
#
# frozen_string_literal: true

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 2,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'set'
require 'test/unit'
require 'bio/dcj'

module Bio
  class TestDCJDistanceAndPath < Test::Unit::TestCase
    def test_reversal_distance_for_three_blocks
      source = Bio::DCJ.genome([[1, 2, 3]])
      target = Bio::DCJ.genome([[1, -2, 3]])

      assert_equal(1, Bio::DCJ.distance(source, target))
      assert_sorts_source_to_target(source, target)
    end

    def test_linear_and_circular_mixture
      source = Bio::DCJ.genome([[1, 2], { blocks: [3, 4], circular: true }])
      target = Bio::DCJ.genome([{ blocks: [2, 3], circular: true }, [4, 1]])

      assert_equal(4, Bio::DCJ.distance(source, target))
      assert_sorts_source_to_target(source, target)
    end

    def test_circular_rotation_and_reversal_do_not_change_distance
      source = Bio::DCJ.genome([{ blocks: [1, 2, 3], circular: true }])
      rotated = Bio::DCJ.genome([{ blocks: [2, 3, 1], circular: true }])
      reversed = Bio::DCJ.genome([{ blocks: [-1, -3, -2], circular: true }])
      target = Bio::DCJ.genome([[1, 2], [-3]])

      assert_equal(source, rotated)
      assert_equal(source, reversed)
      assert_equal(2, Bio::DCJ.distance(source, target))
      assert_equal(2, Bio::DCJ.distance(rotated, target))
      assert_equal(2, Bio::DCJ.distance(reversed, target))
      assert_sorts_source_to_target(rotated, target)
    end

    def test_fusing_two_linear_chromosomes_into_one_circular_chromosome
      source = Bio::DCJ.genome([[1, 2], [3]])
      target = Bio::DCJ.genome([{ blocks: [1, 2, 3], circular: true }])

      assert_equal(2, Bio::DCJ.distance(source, target))
      assert_sorts_source_to_target(source, target)
    end

    def test_single_circular_block_to_linear_block
      source = Bio::DCJ.genome([{ blocks: [1], circular: true }])
      target = Bio::DCJ.genome([[-1]])

      assert_equal(1, Bio::DCJ.distance(source, target))
      assert_sorts_source_to_target(source, target)
    end

    def test_exhaustive_state_graph_for_three_blocks
      source = Bio::DCJ.genome([{ blocks: [1, 2, 3], circular: true }])
      target = Bio::DCJ.genome([[1, 2], [-3]])

      assert_equal(bruteforce_dcj_distance(source, target),
                   Bio::DCJ.distance(source, target))
      assert_sorts_source_to_target(source, target)
    end

    def test_adjacencies_distinguish_block_ends_and_keep_telomeres_half_open
      genome = Bio::DCJ.genome([[1, 2]])

      assert_equal([[[1, :head], [2, :tail]],
                    [[1, :tail]],
                    [[2, :head]]],
                   genome.adjacencies)
    end

    def test_each_returned_operation_must_act_on_current_genome
      source = Bio::DCJ.genome([[1, 2, 3]])
      target = Bio::DCJ.genome([[3, 2, 1]])
      operations = Bio::DCJ.shortest_sorting_operations(source, target)

      current = source
      assert_nothing_raised do
        operations.each do |operation|
          current = Bio::DCJ.apply(current, operation)
        end
      end
      assert_equal(target, current)

      invalid_operation = Bio::DCJ.operation(operations.first.joins,
                                            operations.first.cuts)
      assert_raise(ArgumentError) { source.apply(invalid_operation) }
    end

    def test_missing_duplicate_and_invalid_blocks_are_rejected
      assert_raise(ArgumentError) do
        source = Bio::DCJ.genome([[1, 2, 3]])
        target = Bio::DCJ.genome([[1, 2]])
        Bio::DCJ.distance(source, target)
      end
      assert_raise(ArgumentError) do
        Bio::DCJ.genome([[1, 1, 2]])
      end
      assert_raise(ArgumentError) do
        Bio::DCJ.genome([[0, 2]])
      end
      assert_raise(ArgumentError) do
        Bio::DCJ.genome([['1', 2]])
      end
      assert_raise(ArgumentError) do
        Bio::DCJ.genome([[]])
      end
      assert_raise(ArgumentError) do
        Bio::DCJ.genome([])
      end
      assert_raise(ArgumentError) do
        Bio::DCJ.genome(nil)
      end
      assert_raise(ArgumentError) do
        Bio::DCJ.genome([[1, 2]]).distance_to(Bio::DCJ.genome([[1, 2, 3]]))
      end
    end

    def test_string_block_labels
      source = Bio::DCJ.genome([['+a', '+b']])
      target = Bio::DCJ.genome([['+a', '-b']])

      assert_equal(1, Bio::DCJ.distance(source, target))
      assert_sorts_source_to_target(source, target)
    end

    private

    def assert_sorts_source_to_target(source, target)
      operations = source.shortest_sorting_operations(target)
      assert_equal(Bio::DCJ.distance(source, target), operations.size)

      current = source
      operations.each do |operation|
        assert_kind_of(Bio::DCJ::Operation, operation)
        current = Bio::DCJ.apply(current, operation)
      end
      assert_equal(target, current)
      assert_equal(target.adjacencies, current.adjacencies)
    end

    def bruteforce_dcj_distance(source, target)
      queue = [[source.adjacencies, 0]]
      visited = Set[source.adjacencies]

      until queue.empty?
        state, depth = queue.shift
        return depth if state == target.adjacencies

        neighbour_partitions(state).each do |neighbour|
          next if visited.include?(neighbour)

          visited << neighbour
          queue << [neighbour, depth + 1]
        end
      end

      raise 'brute force search could not reach target'
    end

    def neighbour_partitions(state)
      endpoints = state.flatten(1)
      candidates = []
      internal_edges = state.select { |edge| edge.size == 2 }
      telomere_edges = state.select(&:one?)

      internal_edges.each do |edge|
        candidates << canonical((state - [edge]) + [[edge.first], [edge.last]])
      end

      telomere_edges.combination(2) do |left, right|
        candidates << canonical((state - [left, right]) + [[left.first, right.first]])
      end

      internal_edges.each do |edge|
        telomere_edges.each do |telomere|
          endpoint = telomere.first
          candidates << canonical((state - [edge, telomere]) +
                                  [[edge.first, endpoint], [edge.last]])
          candidates << canonical((state - [edge, telomere]) +
                                  [[edge.first], [edge.last, endpoint]])
        end
      end

      internal_edges.combination(2) do |left, right|
        candidates << canonical((state - [left, right]) +
                                [[left.first, right.first],
                                 [left.last, right.last]])
        candidates << canonical((state - [left, right]) +
                                [[left.first, right.last],
                                 [left.last, right.first]])
      end

      candidates.select { |state| valid_state?(state, endpoints) }.uniq
    end

    def valid_state?(state, endpoints)
      actual = state.flatten(1).sort_by(&:inspect)
      expected = endpoints.sort_by(&:inspect)
      return false unless actual == expected

      state.all? { |edge| edge.size.between?(1, 2) }
    end

    def canonical(state)
      state.map { |edge| edge.sort_by(&:inspect) }.sort_by(&:inspect)
    end
  end
end
