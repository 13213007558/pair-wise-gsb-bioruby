#
# = bio/pcr.rb - sequence-level PCR simulation
#
# License:: The Ruby License
#

# Sequence-level PCR simulation.
module Bio
  require 'bio/sequence'

  # Simulates primer binding and amplicon formation on nucleic acid templates.
  class PCR
    autoload :Coordinates, 'bio/pcr/coordinates'
    autoload :Matcher, 'bio/pcr/matcher'
    autoload :Options, 'bio/pcr/options'
    autoload :Site, 'bio/pcr/result'
    autoload :Product, 'bio/pcr/result'

    def self.simulate(template, forward_primer, reverse_primer, **options)
      new(**options).simulate(template, forward_primer, reverse_primer)
    end

    def initialize(circular: false, max_mismatches: 0,
                   max_total_mismatches: nil,
                   three_prime_mismatches: 0, three_prime_region: 1,
                   min_product_length: 1, max_product_length: nil)
      @circular = circular
      @max_mismatches = max_mismatches
      @max_total_mismatches = max_total_mismatches
      @three_prime_mismatches = three_prime_mismatches
      @three_prime_region = three_prime_region
      @min_product_length = min_product_length
      @max_product_length = max_product_length

      Options.validate(max_mismatches, max_total_mismatches,
                       three_prime_mismatches, three_prime_region,
                       min_product_length, max_product_length)
    end

    def simulate(template, forward_primer, reverse_primer)
      template_sequence = dna(template)
      forward_sequence = dna(forward_primer)
      reverse_sequence = dna(reverse_primer)
      template_length = template_sequence.length
      max_product = @max_product_length || template_length

      Options.validate_sequences(template_sequence, forward_sequence,
                                 reverse_sequence, @circular,
                                 @min_product_length, max_product)

      forward_sites = binding_sites(template_sequence, forward_sequence,
                                    '+')
      reverse_sites = binding_sites(template_sequence, reverse_sequence,
                                    '-')

      site_pairs = forward_sites.product(reverse_sites)
      pairs = site_pairs.filter_map do |forward, reverse|
        length = Coordinates.product_length(template_length,
                                           forward,
                                           reverse,
                                           @circular)
        next unless length.between?(@min_product_length, max_product)

        if @max_total_mismatches &&
           Matcher.total_count(forward, reverse) > @max_total_mismatches
          next
        end

        sequence = Coordinates.product_sequence(template_sequence,
                                               forward,
                                               length)
        Product.new(sequence, forward, reverse)
      end

      pairs.sort_by do |product|
        [product.forward_primer_site.five_prime_position,
         product.reverse_primer_site.five_prime_position,
         product.length]
      end
    end

    private

    def dna(sequence)
      Bio::Sequence::NA.new(sequence.to_s).dna
    end

    def binding_sites(template, primer, strand)
      Matcher.positions(template, primer, @circular).filter_map do |position|
        mismatches = Matcher.mismatches(template, primer, position, strand)
        mismatch_count = mismatches.length
        next unless mismatch_count <= @max_mismatches
        three_prime_count = Matcher.three_prime_count(mismatches, primer,
                                                       @three_prime_region)
        next unless three_prime_count <= @three_prime_mismatches

        Coordinates.site(position, primer, strand, mismatches,
                         template.length, @circular)
      end
    end
  end
end
