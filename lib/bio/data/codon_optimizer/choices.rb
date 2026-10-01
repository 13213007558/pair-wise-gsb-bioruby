# frozen_string_literal: true

require 'bio/data/codon_optimizer/support'

module Bio
  class CodonOptimizer
    # Builds and validates synonymous codon choices for one coding sequence.
    class Choices

      def initialize(codon_table, sequence, fixed_positions:,
                     fix_start:, fix_stop:)
        @codon_table = codon_table
        @rna = sequence.downcase.include?('u')
        @sequence = sequence.downcase.tr('u', 't')
        @fixed_positions = fixed_positions
        @fix_start = fix_start
        @fix_stop = fix_stop
      end

      def conflicts
        validator = SequenceValidator.new(@sequence, @fixed_positions)
        @conflicts ||= validator.conflicts
      end

      def lists
        return [] unless conflicts.empty?

        @builder ||= ChoiceBuilder.new(
          @codon_table,
          codons,
          @fixed_positions,
          @fix_start,
          @fix_stop
        )
        @builder.build
      end

      def build_conflicts
        lists
        @builder&.conflicts || []
      end

      def rna?
        @rna
      end

      private

      def codons
        @sequence.scan(/.{3}/)
      end

    end

    class SequenceValidator

      def initialize(sequence, fixed_positions)
        @sequence = sequence
        @fixed_positions = fixed_positions
      end

      def conflicts
        sequence_conflicts + position_conflicts
      end

      private

      def sequence_conflicts
        result = []
        unless @sequence.length.remainder(3).zero?
          result << Support.conflict(
            'invalid_sequence_length',
            'coding sequence length must be a positive multiple of 3'
          )
        end
        return result if @sequence.match?(/\A[atgc]+\z/)

        result << Support.conflict(
          'invalid_nucleotide',
          'coding sequence must contain only A, C, G, T, or U'
        )
        result
      end

      def position_conflicts
        count = @sequence.length / 3
        out_of_range = @fixed_positions.reject { |position| position < count }
        out_of_range.map do |position|
          Support.conflict(
            'fixed_position_out_of_range',
            "fixed position #{position} is outside the coding sequence",
            position: position
          )
        end
      end

    end

    class ChoiceBuilder

      attr_reader :conflicts

      def initialize(codon_table, codons, fixed_positions,
                     fix_start, fix_stop)
        @codon_table = codon_table
        @codons = codons
        @fixed_positions = fixed_positions
        @fix_start = fix_start
        @fix_stop = fix_stop
        @conflicts = []
      end

      def build
        lists = @codons.each_with_index.map do |codon, position|
          choices_for(codon, position)
        end
        @conflicts.empty? ? lists : []
      end

      private

      def choices_for(codon, position)
        amino_acid = @codon_table.table[codon]
        return invalid_codon(codon, position) if amino_acid.nil?

        synonyms = synonyms_for(amino_acid)
        if synonyms.empty?
          return missing_amino_acid(amino_acid, position)
        end
        return synonyms unless fixed_positions.include?(position)
        return [codon] if synonyms.include?(codon)

        @conflicts << fixed_codon_conflict(codon, amino_acid, position)
        []
      end

      def fixed_positions
        fixed = @fixed_positions.select { |position| position < @codons.length }
        fixed << 0 if @fix_start
        if @fix_stop && @codon_table.table[@codons.last] == '*'
          fixed << @codons.length - 1
        end
        fixed
      end

      def synonyms_for(amino_acid)
        @codon_table.table.each_key.select do |codon|
          codon.match?(/\A[atgc]{3}\z/) &&
            @codon_table.table[codon] == amino_acid
        end.sort
      end

      def invalid_codon(codon, position)
        @conflicts << Support.conflict(
          'invalid_codon',
          "#{codon} at position #{position} is not defined by the codon table",
          position: position,
          codon: codon
        )
        []
      end

      def missing_amino_acid(amino_acid, position)
        @conflicts << Support.conflict(
          'unsupported_amino_acid',
          "no codon encodes #{amino_acid} in the supplied codon table",
          position: position
        )
        []
      end

      def fixed_codon_conflict(codon, amino_acid, position)
        Support.conflict(
          'fixed_codon_mismatch',
          "fixed codon #{codon} does not encode #{amino_acid}",
          position: position,
          codon: codon
        )
      end

    end

  end

end
