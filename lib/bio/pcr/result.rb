# frozen_string_literal: true

# Result objects returned by Bio::PCR.
module Bio
  class PCR
    # A primer binding site accepted by the simulator.
    Site = Struct.new(:primer, :strand, :five_prime_position,
                      :three_prime_position, :mismatch_positions,
                      :mismatches)

    # An amplicon generated from an oriented pair of binding sites.
    Product = Struct.new(:sequence, :forward_primer_site,
                         :reverse_primer_site) do
      def length
        sequence.length
      end

      def mismatch_positions
        {
          forward_primer: forward_primer_site.mismatch_positions,
          reverse_primer: reverse_primer_site.mismatch_positions
        }
      end

      def mismatches
        {
          forward_primer: forward_primer_site.mismatches,
          reverse_primer: reverse_primer_site.mismatches
        }
      end
    end
  end
end
