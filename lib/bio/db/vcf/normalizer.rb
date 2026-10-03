# frozen_string_literal: true

module Bio
  class VCF
    # Left-aligns and minimizes explicit sequence alleles against a reference.
    class Normalizer
      def self.symbolic_allele?(allele)
        allele.start_with?('<') || allele == '*' ||
          allele.include?('[') || allele.include?(']')
      end

      def initialize(reference = nil, sequence_start: 1)
        @reference = reference
        @sequence_start = sequence_start
      end

      def normalize(record, split: false)
        if record.symbolic?
          mark_unprocessed(record.dup)
        elsif split && record.alt.length > 1
          record.alt.each_index.collect do |index|
            normalize_explicit(split_record(record, index))
          end
        else
          normalize_explicit(record.dup)
        end
      end

      alias normalize_record normalize

      def split(record)
        result = record.dup
        return [mark_unprocessed(result)] if result.symbolic?
        return [result] unless result.alt.length > 1

        result.alt.each_index.collect { |index| split_record(result, index) }
      end

      private

      def mark_unprocessed(record)
        record.normalization_status = :unprocessed
        record
      end

      def normalize_explicit(record)
        return mark_unprocessed(record) if record.symbolic?

        sequence = reference_sequence(record.chrom)
        unless sequence
          raise ArgumentError, "missing reference sequence for #{record.chrom}"
        end
        apply_reference(record, sequence)
      end

      def apply_reference(record, sequence)
        AlleleNormalizer.new(sequence, @sequence_start).normalize(record)
      end

      def split_record(record, selected_index)
        RecordSplitter.new.split(record, selected_index)
      end

      def reference_sequence(chrom)
        reference = @reference.is_a?(Hash) ? @reference[chrom] : @reference
        reference && reference.to_s
      end
    end
  end
end
