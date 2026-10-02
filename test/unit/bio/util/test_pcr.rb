# frozen_string_literal: true
#
# test/unit/bio/util/test_pcr.rb - Unit test for Bio::PCR.
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/util/pcr'

module Bio
  class TestPCR < Test::Unit::TestCase
    # Forward and reverse primers used throughout the tests.
    # Neither is a reverse-complement palindrome, so the strand of
    # every binding site is unambiguous.
    FORWARD = 'ACGTTGCA' # reverse complement: TGCAACGT
    REVERSE = 'GGCATTGC' # reverse complement: GCAATGCC

    # Template with a single perfect product:
    #   FORWARD site at 5...13, reverse-complement site of REVERSE
    #   at 25...33.
    TEMPLATE = 'CCCCC' + 'ACGTTGCA' + 'TTAGTTAGTTAG' + 'GCAATGCC' + 'GGGGG'
    PRODUCT_SEQ = 'ACGTTGCA' + 'TTAGTTAGTTAG' + 'GCAATGCC'

    def test_new
      assert(Bio::PCR.new(TEMPLATE))
      assert(Bio::PCR.new(Bio::Sequence::NA.new(TEMPLATE)))
      assert(Bio::PCR.new(TEMPLATE, circular: true))
    end

    def test_new_invalid
      assert_raise(ArgumentError) { Bio::PCR.new('') }
    end

    def test_unique_perfect_product
      products = Bio::PCR.new(TEMPLATE).amplify(FORWARD, REVERSE)
      assert_equal(1, products.size)

      product = products.first
      assert_equal(5, product.start)
      assert_equal(33, product.stop)
      assert_equal(28, product.length)
      assert_equal(:plus, product.strand)
      assert_kind_of(Bio::Sequence::NA, product.sequence)
      assert_equal(PRODUCT_SEQ, product.sequence.to_s.upcase)

      assert_equal(5, product.forward_match.start)
      assert_equal(:plus, product.forward_match.strand)
      assert_equal([], product.forward_match.mismatches)
      assert_equal(25, product.reverse_match.start)
      assert_equal(:minus, product.reverse_match.strand)
      assert_equal([], product.reverse_match.mismatches)
    end

    def test_false_positive_site_with_only_three_prime_error
      # A second reverse-primer site differing only at the 3'-terminal
      # base of the primer (the first base of its reverse complement).
      template = TEMPLATE[0...33] + 'TTAG' + 'ACAATGCC' + 'TT'
      pcr = Bio::PCR.new(template)

      # Total mismatches are tolerated, but the strict 3'-end rule
      # rejects the false-positive site.
      products = pcr.amplify(FORWARD, REVERSE, max_mismatches: 2)
      assert_equal(1, products.size)
      assert_equal(33, products.first.stop)

      # Relaxing the 3'-end rule reveals the false-positive product.
      products = pcr.amplify(FORWARD, REVERSE,
                             max_mismatches: 2,
                             max_three_prime_mismatches: 1)
      assert_equal(2, products.size)
      false_positive = products.last
      assert_equal(45, false_positive.stop)
      assert_equal([7], false_positive.reverse_match.mismatches)
      assert_equal(1, false_positive.reverse_match.three_prime_mismatches(5))
    end

    def test_degenerate_iupac_primer
      template = 'ACGTTGCA' + 'TTAG' + 'ACGTTGCG' + 'TTAG' + 'GCAATGCC'
      products = Bio::PCR.new(template).amplify('ACGTTGCR', REVERSE)
      assert_equal([32, 20], products.map(&:length))
      assert_equal([0, 12], products.map(&:start))
    end

    def test_mismatch_positions_of_forward_primer
      template = 'CC' + 'AGGTTGCA' + 'TTAG' + 'GCAATGCC' + 'CC'
      pcr = Bio::PCR.new(template)

      assert_equal([], pcr.amplify(FORWARD, REVERSE))

      products = pcr.amplify(FORWARD, REVERSE, max_mismatches: 1)
      assert_equal(1, products.size)
      assert_equal([1], products.first.forward_match.mismatches)
      assert_equal(0, products.first.forward_match.three_prime_mismatches(5))
    end

    def test_three_prime_rule_on_forward_primer
      # Single mismatch at the 3'-terminal base of the forward primer.
      template = 'CC' + 'ACGTTGCG' + 'TTAG' + 'GCAATGCC' + 'CC'
      pcr = Bio::PCR.new(template)

      assert_equal([], pcr.amplify(FORWARD, REVERSE, max_mismatches: 1))

      products = pcr.amplify(FORWARD, REVERSE,
                             max_mismatches: 1,
                             max_three_prime_mismatches: 1)
      assert_equal(1, products.size)
      assert_equal([7], products.first.forward_match.mismatches)
    end

    def test_product_length_range
      template = 'ACGTTGCA' + 'TTAG' + 'ACGTTGCG' + 'TTAG' + 'GCAATGCC'
      pcr = Bio::PCR.new(template)

      products = pcr.amplify('ACGTTGCR', REVERSE, min_product_length: 25)
      assert_equal([32], products.map(&:length))

      products = pcr.amplify('ACGTTGCR', REVERSE, max_product_length: 25)
      assert_equal([20], products.map(&:length))

      products = pcr.amplify('ACGTTGCR', REVERSE,
                             min_product_length: 21,
                             max_product_length: 31)
      assert_equal([], products)
    end

    def test_product_on_minus_strand
      # REVERSE matches the plus strand directly and FORWARD matches
      # the minus strand downstream of it.
      template = 'CC' + 'GGCATTGC' + 'TTAGTTAG' + 'TGCAACGT' + 'GG'
      products = Bio::PCR.new(template).amplify(FORWARD, REVERSE)
      assert_equal(1, products.size)

      product = products.first
      assert_equal(2, product.start)
      assert_equal(26, product.stop)
      assert_equal(:minus, product.strand)
      assert_equal('ACGTTGCA' + 'CTAACTAA' + 'GCAATGCC',
                   product.sequence.to_s.upcase)
      assert_equal(:minus, product.forward_match.strand)
      assert_equal(18, product.forward_match.start)
      assert_equal(:plus, product.reverse_match.strand)
      assert_equal(2, product.reverse_match.start)
    end

    def test_divergent_sites_yield_no_product
      # Only left-facing sites (primers would extend away from
      # each other).
      template = 'TGCAACGT' + 'TTAG' + 'GCAATGCC'
      assert_equal([], Bio::PCR.new(template).amplify(FORWARD, REVERSE))

      # Only right-facing sites.
      template = 'ACGTTGCA' + 'TTAG' + 'GGCATTGC'
      assert_equal([], Bio::PCR.new(template).amplify(FORWARD, REVERSE))
    end

    def test_circular_template_wraps_at_most_once
      # On a linear template the sites face away from each other;
      # on a circular template the product wraps around the origin.
      template = 'GCAATGCC' + 'TTAGTTAG' + 'ACGTTGCA'
      assert_equal([], Bio::PCR.new(template).amplify(FORWARD, REVERSE))

      products = Bio::PCR.new(template, circular: true).amplify(FORWARD, REVERSE)
      assert_equal(1, products.size)

      product = products.first
      assert_equal(16, product.start)
      assert_equal(32, product.stop)
      assert_equal(16, product.length)
      assert_operator(product.length, :<=, template.length)
      assert_equal('ACGTTGCA' + 'GCAATGCC', product.sequence.to_s.upcase)
    end

    def test_invalid_primer
      pcr = Bio::PCR.new(TEMPLATE)
      assert_raise(ArgumentError) { pcr.amplify('', REVERSE) }
      assert_raise(ArgumentError) { pcr.amplify('ACGTXGCA', REVERSE) }
    end

    def test_unknown_option
      pcr = Bio::PCR.new(TEMPLATE)
      assert_raise(ArgumentError) { pcr.amplify(FORWARD, REVERSE, mismatch: 1) }
    end
  end
end
