#
# = test/unit/bio/test_pcr.rb - Unit test for Bio::PCR
#
# License:: The Ruby License
#

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 2,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio'

module Bio
  class TestPCR < Test::Unit::TestCase
    def setup
      @reverse_primer = 'accatgctaa'
      @template = 'aaa' \
                  'atgcgatacg' \
                  'atgcgatact' \
                  'ttagcatggt' \
                  'ccc'
    end

    def test_unique_exact_product
      products = Bio::PCR.simulate(@template, 'atgcgatact',
                                   @reverse_primer,
                                   min_product_length: 10,
                                   max_product_length: 40)

      assert_equal(1, products.length)
      assert_equal('atgcgatactttagcatggt', products.first.sequence)
      assert_equal(20, products.first.length)
      assert_equal('+', products.first.forward_primer_site.strand)
      assert_equal('-', products.first.reverse_primer_site.strand)
      assert_equal(14, products.first.forward_primer_site.five_prime_position)
      assert_equal(33, products.first.reverse_primer_site.five_prime_position)
      assert_equal([], products.first.mismatch_positions[:forward_primer])
      assert_equal([], products.first.mismatch_positions[:reverse_primer])
    end

    def test_three_prime_terminal_false_positive_is_rejected
      strict = Bio::PCR.new(min_product_length: 10,
                           max_product_length: 40)
      permissive = Bio::PCR.new(max_mismatches: 1,
                                three_prime_mismatches: 1,
                                min_product_length: 10,
                                max_product_length: 40)

      assert_equal(1, strict.simulate(@template, 'atgcgatact',
                                      @reverse_primer).length)
      products = permissive.simulate(@template, 'atgcgatact',
                                     @reverse_primer)

      assert_equal([30, 20], products.map(&:length))
      mismatching_product = products.find { |product| product.length == 30 }
      assert_equal([10],
                   mismatching_product.mismatch_positions[:forward_primer])
    end

    def test_internal_mismatch_positions
      products = Bio::PCR.simulate('atgcgatacgttagcatggt',
                                   'atgcaatacg', 'accatggtaa',
                                   max_mismatches: 1,
                                   min_product_length: 10,
                                   max_product_length: 30)

      assert_equal(1, products.length)
      assert_equal([5], products.first.mismatch_positions[:forward_primer])
      assert_equal([7], products.first.mismatch_positions[:reverse_primer])
    end

    def test_total_mismatch_limit_spans_both_primers
      template = 'atgcaagacgttagcatggt'
      options = {
        max_mismatches: 1,
        three_prime_mismatches: 0,
        min_product_length: 10,
        max_product_length: 30
      }

      assert_empty(Bio::PCR.simulate(template, 'atgcaatacg',
                                    'accatggtaa',
                                    max_total_mismatches: 1,
                                    **options))
      assert_equal(1, Bio::PCR.simulate(template,
                                        'atgcaatacg',
                                        'accatggtaa',
                                        max_total_mismatches: 2,
                                        **options).length)
    end

    def test_iupac_degenerate_primers_match
      products = Bio::PCR.simulate('atgcgttagc', 'atgcr', 'gctar')

      assert_equal(1, products.length)
      assert_equal('atgcgttagc', products.first.sequence)
      assert_equal([], products.first.mismatch_positions[:forward_primer])
      assert_equal([], products.first.mismatch_positions[:reverse_primer])
    end

    def test_primers_must_face_inward
      template = 'ttagcatggtatgcgatacg'

      assert_empty(Bio::PCR.simulate(template, 'atgcgatacg',
                                     @reverse_primer))
    end

    def test_repeated_binding_sites_keep_distinct_products
      template = 'atgcNNNNNtaacNNNNNatgcNNNNNtaac'
      products = Bio::PCR.simulate(template, 'atgc', 'gtta')

      assert_equal([13, 31, 13], products.map(&:length))
      assert_equal([1, 1, 19],
                   products.map { |product|
                     product.forward_primer_site.five_prime_position
                   })
      assert_equal([13, 31, 31],
                   products.map { |product|
                     product.reverse_primer_site.five_prime_position
                   })
    end

    def test_product_length_range
      products = Bio::PCR.simulate(@template, 'atgcgatact',
                                   @reverse_primer,
                                   min_product_length: 25,
                                   max_product_length: 35,
                                   max_mismatches: 1,
                                   three_prime_mismatches: 1)

      assert_equal([30], products.map(&:length))
    end

    def test_circular_template_can_wrap_at_most_once
      template = 'gcat' + ('c' * 22) + 'atgc'
      linear = Bio::PCR.simulate(template, 'atgc', 'atgc')
      circular = Bio::PCR.simulate(template, 'atgc', 'atgc',
                                   circular: true)

      assert_empty(linear)
      assert_equal(1, circular.length)
      assert_equal(8, circular.first.length)
      assert_equal('atgcgcat', circular.first.sequence)
      assert_equal(27, circular.first.forward_primer_site.five_prime_position)
      assert_equal(4, circular.first.reverse_primer_site.five_prime_position)
    end

    def test_circular_full_length_product
      template = 'a' * 10
      products = Bio::PCR.simulate(template, 'aaaaa', 'ttttt',
                                   circular: true)

      assert_equal(10, products.length)
      assert(products.all? { |product| product.length == 10 })
    end

    def test_sequence_nas_have_a_pcr_convenience_method
      sequence = Bio::Sequence::NA.new(@template)
      products = sequence.pcr('atgcgatact', @reverse_primer,
                              min_product_length: 10,
                              max_product_length: 40)

      assert_equal(1, products.length)
      assert_equal('atgcgatactttagcatggt', products.first.sequence)
    end
  end
end
