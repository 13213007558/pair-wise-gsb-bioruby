# frozen_string_literal: true

# bio/codon_optimizer.rb - exact synonymous codon optimizer

require 'bio/data/codontable'
require 'bio/sequence'
require 'bio/sequence/na'
require_relative 'codon_optimizer/validator'
require_relative 'codon_optimizer/automaton'
require_relative 'codon_optimizer/problem'
require_relative 'codon_optimizer/solver'

module Bio
  # Exact optimizer for synonymous coding sequences.
  #
  # An Aho--Corasick automaton detects forbidden substrings, including
  # matches that cross codon boundaries.  A dynamic program keeps every
  # reachable automaton state and GC count, instead of greedily choosing one
  # codon first and attempting random repairs afterward.
  #
  # Positions in +fixed_positions+ are zero-based codon positions.  Values in
  # +codon_weights+ are per-codon scores; the optimizer maximizes their sum,
  # with larger values preferred.  By default, a start codon at
  # position 0 and a stop codon at the final position are preserved because
  # initiation codons such as GTG do not encode their initiating amino acid
  # under normal translation.  Pass +start_codon_strategy: :optimize+ or
  # +stop_codon_strategy: :optimize+ to allow the corresponding position to
  # choose a normal synonymous codon from the caller-supplied table.
  class CodonOptimizer
    # Optimization result returned by #optimize.
    Result = Struct.new(:sequence, :score, :gc_count, :conflicts,
                        keyword_init: true) do
      def feasible?
        !sequence.nil?
      end
    end

    private_constant :Automaton, :Problem, :Solver, :Validator

    # Returns the highest-scoring sequence satisfying all constraints.
    def self.optimize(sequence, table, **options)
      new(table, **options).optimize(sequence)
    end

    def initialize(table, **options)
      @problem = Problem.new(table, options)
      @solver = Solver.new(@problem)
    end

    def optimize(sequence)
      source = Validator.normalize_nucleotides(sequence, 'sequence')
      raise ArgumentError, 'sequence length must be a multiple of 3' \
        unless (source.length % 3).zero?

      codons = source.scan(/.{3}/)
      @problem.validate_sequence(codons)
      automaton = Automaton.new(@problem.patterns)
      layer = @solver.run(codons, automaton)
      best = @solver.best(layer)
      return infeasible_result(codons, automaton) unless best

      build_result(sequence, best)
    end

    private

    def infeasible_result(codons, automaton)
      forbidden_layer = @solver.run(codons, automaton, nil)
      conflicts = conflicts(forbidden_layer)
      Result.new(sequence: nil, score: nil, gc_count: nil,
                 conflicts: conflicts)
    end

    def conflicts(layer)
      return [{ 'type' => 'forbidden_substring' }] if layer.empty?

      achievable = layer.values.map(&:gc_count).minmax
      [{
        'type' => 'gc_count',
        'requested_min' => @problem.min_gc_count,
        'requested_max' => @problem.max_gc_count,
        'achievable_min' => achievable[0],
        'achievable_max' => achievable[1]
      }]
    end

    def build_result(source_input, node)
      optimized = @solver.reconstruct(node)
      output = convert_output(source_input, optimized)
      Result.new(sequence: output, score: node.score,
                 gc_count: node.gc_count, conflicts: [])
    end

    def convert_output(source_input, dna_sequence)
      output = if source_input.to_str.downcase.include?('u')
                 dna_sequence.dup.tr('t', 'u')
               else
                 dna_sequence
               end
      Bio::Sequence::NA.new(output)
    end
  end
end
