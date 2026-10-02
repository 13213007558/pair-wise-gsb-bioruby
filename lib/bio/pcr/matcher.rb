# frozen_string_literal: true

# IUPAC-aware matching helpers used by Bio::PCR.
module Bio
  class PCR
    # IUPAC-aware primer matching on either strand.
    module Matcher
      module_function

      def positions(template, primer, circular)
        return 0...template.length if circular
        return [] if primer.length > template.length

        0..(template.length - primer.length)
      end

      def mismatches(template, primer, position, strand)
        aligned_primer = strand == '+' ? primer : primer.reverse
        aligned_primer.each_char.with_index.filter_map do |base, offset|
          template_base = template[(position + offset) % template.length]
          target_base = target_base(template_base, strand)
          next if compatible?(base, target_base)

          strand == '+' ? offset + 1 : primer.length - offset
        end
      end

      def three_prime_count(mismatches, primer, region)
        first_position = primer.length - region + 1
        mismatches.count { |position| position >= first_position }
      end

      def total_count(forward, reverse)
        forward.mismatches + reverse.mismatches
      end

      def target_base(base, strand)
        return base if strand == '+'

        Bio::Sequence::NA.new(base).forward_complement
      end

      def compatible?(primer_base, template_base)
        Bio::NucleicAcid.to_re(primer_base).match?(template_base)
      end
    end
  end
end
