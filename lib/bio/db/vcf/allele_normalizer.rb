# frozen_string_literal: true

module Bio
  class VCF
    # Normalizes explicit sequence alleles against one reference sequence.
    class AlleleNormalizer
      def initialize(sequence, sequence_start)
        @sequence = sequence
        @sequence_start = sequence_start
      end

      def normalize(record)
        alleles = [record.ref] + record.alt
        alleles.collect!(&:upcase)
        offset = record.pos - @sequence_start
        validate_reference(record, alleles, offset)

        offset = left_align(alleles, offset)
        record.pos = offset + @sequence_start
        record.ref = alleles.first
        record.alt = alleles[1..-1]
        record.normalization_status = :normalized
        record
      end

      private

      def validate_reference(record, alleles, offset)
        expected = @sequence[offset, alleles.first.length]
        return if expected && expected.upcase == alleles.first

        raise ArgumentError, 'REF does not match the reference sequence'
      end

      def left_align(alleles, offset)
        offset = trim_prefix(alleles, offset)
        loop do
          trim_suffix(alleles, offset)
          break unless offset > 0 && alleles.any?(&:empty?)

          prepend_previous_base(alleles, offset)
          offset -= 1
          offset = trim_prefix(alleles, offset)
        end
        offset
      end

      def prepend_previous_base(alleles, offset)
        previous_base = @sequence[offset - 1].upcase
        alleles.each { |allele| allele.replace(previous_base + allele) }
      end

      def trim_suffix(alleles, offset)
        while alleles.none?(&:empty?) && shared_suffix?(alleles) &&
              suffix_trim_allowed?(alleles, offset)
          alleles.each(&:chop!)
        end
      end

      def shared_suffix?(alleles)
        alleles.collect { |allele| allele[-1] }.uniq.length == 1
      end

      def suffix_trim_allowed?(alleles, offset)
        offset > 0 || alleles.all? { |allele| allele.length > 1 }
      end

      def trim_prefix(alleles, offset)
        while alleles.all? { |allele| allele.length > 1 } && shared_prefix?(alleles)
          alleles.each { |allele| allele.slice!(0) }
          offset += 1
        end
        offset
      end

      def shared_prefix?(alleles)
        alleles.collect { |allele| allele[0] }.uniq.length == 1
      end
    end
  end
end
