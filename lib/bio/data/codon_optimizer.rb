# frozen_string_literal: true
#
# = bio/data/codon_optimizer.rb - Exact synonymous codon optimizer
#
# Copyright:: Copyright (C) 2026 BioRuby project
# License:: The Ruby License
#

require 'bio/data/codontable'
require 'bio/data/codon_optimizer/support'
require 'bio/data/codon_optimizer/configuration'
require 'bio/data/codon_optimizer/choices'
require 'bio/data/codon_optimizer/search'

module Bio
  # Exact optimizer for synonymous codon choices.
  #
  # Search and start/stop policy details are documented on Search and
  # Configuration.
  class CodonOptimizer
    # Result returned by CodonOptimizer#optimize.
    Result = Struct.new(:sequence, :score, :gc_count, :conflicts,
                        keyword_init: true) do

      def feasible?
        conflicts.empty?
      end
    end

    # The codon table may be numeric or a Bio::CodonTable.
    def initialize(codon_table:, codon_weights: nil,
                   forbidden_substrings: [], fixed_positions: [],
                   gc_range: nil, fix_start: true, fix_stop: true)
      @configuration = Configuration.new(
        codon_table: codon_table,
        codon_weights: Arguments.weights(codon_weights),
        forbidden_substrings: Arguments.patterns(forbidden_substrings),
        fixed_positions: Arguments.positions(fixed_positions),
        gc_range: Arguments.gc_range(gc_range),
        fix_start: fix_start,
        fix_stop: fix_stop
      )
    end

    # Returns an optimized coding sequence, or a Result containing conflicts
    # when no assignment satisfies every constraint.
    def optimize(sequence)
      unless sequence.is_a?(String)
        conflict = Support.conflict('invalid_sequence',
                                    'coding sequence must be a String')
        return Support.infeasible(Result, [conflict])
      end
      return infeasible(empty_sequence_conflict) if sequence.empty?

      choices = @configuration.choices_for(sequence)
      return infeasible(choices.conflicts) unless choices.conflicts.empty?
      lists = choices.lists
      build_conflicts = choices.build_conflicts
      return infeasible(build_conflicts) unless build_conflicts.empty?

      gc_conflict = @configuration.gc_conflict(lists)
      return infeasible([gc_conflict]) if gc_conflict

      result, program_conflict = @configuration.search(lists).run
      return infeasible([program_conflict]) unless result

      Support.feasible(Result, result, choices.rna?)
    end

    # Convenience constructor and optimizer.
    def self.optimize(sequence, codon_table:, **options)
      new(codon_table: codon_table, **options).optimize(sequence)
    end

    private

    def infeasible(conflicts)
      Support.infeasible(Result, conflicts)
    end

    def empty_sequence_conflict
      Support.conflict(
        'invalid_sequence_length',
        'coding sequence length must be a positive multiple of 3'
      )
    end

  end

end
