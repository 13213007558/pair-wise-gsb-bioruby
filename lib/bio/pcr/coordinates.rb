# frozen_string_literal: true

# Coordinate helpers used by Bio::PCR.
module Bio
  class PCR
    # Coordinate conversion for primer sites and products.
    module Coordinates
      module_function

      def site(position, primer, strand, mismatches, template_length,
               circular)
        three_prime_offset = primer.length - 1

        if strand == '+'
          five_prime = position + 1
          three_prime = terminal_position(position, three_prime_offset,
                                          template_length, circular,
                                          primer.length)
        else
          five_prime = terminal_position(position, three_prime_offset,
                                         template_length, circular,
                                         primer.length)
          three_prime = position + 1
        end

        Site.new(primer, strand, five_prime, three_prime, mismatches,
                 mismatches.length)
      end

      def product_length(template_length, forward, reverse, circular)
        forward_start = forward.five_prime_position - 1
        reverse_start = reverse.three_prime_position - 1
        if circular && reverse_start <= forward_start
          reverse_start += template_length
        end

        length = reverse_start + reverse.primer.length - forward_start
        minimum_length = forward.primer.length + reverse.primer.length
        length >= minimum_length ? length : 0
      end

      def product_sequence(template, forward, length)
        start = forward.five_prime_position - 1
        sequence = Array.new(length) do |offset|
          template[(start + offset) % template.length]
        end
        Bio::Sequence::NA.new(sequence.join)
      end

      def terminal_position(position, offset, template_length, circular,
                            primer_length)
        return position + primer_length unless circular

        ((position + offset) % template_length) + 1
      end
    end
  end
end
