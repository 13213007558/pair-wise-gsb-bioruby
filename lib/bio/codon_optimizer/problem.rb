# frozen_string_literal: true

# bio/codon_optimizer/problem.rb - validated optimizer input

module Bio
  class CodonOptimizer
    # Validated caller-supplied optimizer options.
    class Problem
      attr_reader :table, :patterns, :weights, :fixed_positions,
                  :min_gc_count, :max_gc_count, :start_codon_strategy,
                  :stop_codon_strategy

      def initialize(table, options)
        source_table = table.is_a?(CodonTable) ? table : CodonTable[table]
        @table = Marshal.load(Marshal.dump(source_table))
        @patterns = Validator.normalize_patterns(
          options.fetch(:forbidden_substrings, [])
        )
        @weights = Validator.normalize_weights(options.fetch(:codon_weights, {}))
        @fixed_positions = Validator.normalize_positions(
          options.fetch(:fixed_positions, nil)
        )
        @min_gc_count, @max_gc_count = Validator.normalize_gc_range(
          options.fetch(:gc_count, options[:gc_range]),
          options.fetch(:min_gc_count, nil),
          options.fetch(:max_gc_count, nil)
        )
        @start_codon_strategy = Validator.validate_strategy(
          options.fetch(:start_codon_strategy, :preserve),
          :start_codon_strategy
        )
        @stop_codon_strategy = Validator.validate_strategy(
          options.fetch(:stop_codon_strategy, :preserve),
          :stop_codon_strategy
        )
      end

      def validate_sequence(codons)
        Validator.validate_codons(table, codons)
        validate_gc_bounds(codons.length * 3)
        validate_fixed_positions(codons.length)
      end

      def fixed?(codon, index, codon_count)
        fixed_positions.include?(index) ||
          fixed_start?(codon, index) ||
          fixed_stop?(codon, index, codon_count)
      end

      private

      def validate_gc_bounds(maximum_possible)
        if min_gc_count > maximum_possible
          raise ArgumentError, 'min_gc_count exceeds sequence length'
        end
        return unless max_gc_count && max_gc_count > maximum_possible

        raise ArgumentError, 'max_gc_count exceeds sequence length'
      end

      def validate_fixed_positions(codon_count)
        fixed_positions.each do |position|
          next if position < codon_count

          raise ArgumentError, "fixed codon position #{position} is outside " \
                               "the sequence (0...#{codon_count})"
        end
      end

      def fixed_start?(codon, index)
        index.zero? && start_codon_strategy == :preserve &&
          table.start_codon?(codon)
      end

      def fixed_stop?(codon, index, codon_count)
        index == codon_count - 1 && stop_codon_strategy == :preserve &&
          table.stop_codon?(codon)
      end
    end
  end
end
