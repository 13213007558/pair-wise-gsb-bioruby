# frozen_string_literal: true

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 4,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/util/restriction_enzyme'

module Bio
  module TestRestrictionEnzyme
    class TestDigest < Test::Unit::TestCase
      ECO_RI = {
        recognition: 'gaattc',
        cuts: [[0, 4]],
        sensitive_to: [{strand: :+, position: 0, modification: 'm6a'}]
      }

      def test_linear_sticky_digest_reports_cuts_fragments_and_ends
        result = Bio::RestrictionEnzyme.digest('ggaattc', ECO_RI)

        assert_equal(1, result.cuts.size)
        cut = result.cuts.first
        assert_equal(1, cut.primary_cut)
        assert_equal(5, cut.complement_cut)
        assert_equal(:five_prime, cut.sticky_end.kind)
        assert_equal(4, cut.sticky_end.length)
        assert_equal('aatt', cut.sticky_end.sequence)
        assert_equal(['gg    ', '    aattc'], result.fragments.map(&:primary))
        complements = result.fragments.map(&:complement)
        assert_equal(['ccttaa', 'g        '], complements)
        assert(result.fragments.all?(&:linear?))
      end

      def test_circular_blunt_single_cut_is_one_linear_fragment
        rule = {recognition: 'gaattc', cuts: [[2, 2]]}
        result = Bio::RestrictionEnzyme.digest('gaattc', rule,
                                               circular: true)

        assert_equal([2], result.cuts.map(&:primary_cut))
        assert_equal(1, result.fragments.size)
        fragment = result.fragments.first
        assert(fragment.linear?)
        assert_equal('ttcgaa', fragment.primary)
        assert_equal('aagctt', fragment.complement)
      end

      def test_recognition_site_can_cross_the_first_base_of_circular_dna
        result = Bio::RestrictionEnzyme.digest('attcga', ECO_RI,
                                              circular: true)

        assert_equal([4], result.cuts.map(&:primary_cut))
        assert_equal([2], result.cuts.map(&:complement_cut))
        assert_equal(1, result.fragments.size)
        assert(fragment = result.fragments.first)
        assert(fragment.linear?)
        assert_equal('aattcg', fragment.primary.tr(' ', ''))
        assert_equal('gcttaa', fragment.complement.tr(' ', ''))
      end

      def test_rotation_is_equivalent_after_coordinate_mapping
        original = Bio::RestrictionEnzyme.digest('attcga', ECO_RI,
                                                circular: true)
        rotated_sequence = rotate('attcga', 2)
        rotated = Bio::RestrictionEnzyme.digest(rotated_sequence, ECO_RI,
                                               circular: true)
        size = 6

        original_primary = original.fragments.first.primary.tr(' ', '')
        rotated_primary = rotated.fragments.first.primary.tr(' ', '')
        assert_equal(original_primary, rotated_primary)
        assert_equal(original.cuts.map { |cut| (cut.primary_cut - 2) % size },
                     rotated.cuts.map(&:primary_cut))
        mapped_complement = original.cuts.map do |cut|
          (cut.complement_cut - 2) % size
        end
        assert_equal(mapped_complement, rotated.cuts.map(&:complement_cut))
      end

      def test_strand_specific_methylation_blocks_cross_origin_site
        sequence = 'attcga'
        blocked = Bio::RestrictionEnzyme.digest(sequence,
                                               ECO_RI,
                                               circular: true,
                                                methylation: [{
                                                  strand: :+,
                                                  position: 4,
                                                  modification: 'm6a'
                                                }])
        allowed = Bio::RestrictionEnzyme.digest(sequence,
                                               ECO_RI,
                                               circular: true,
                                                methylation: [{
                                                  strand: :-,
                                                  position: 4,
                                                  modification: 'm6a'
                                                }])

        assert(blocked.cuts.empty?)
        assert_equal(1, blocked.blocked_sites.size)
        site = blocked.blocked_sites.first
        assert_equal(4, site.start)
        assert_equal(3, site.end)
        assert_equal([[:+, 4, 'm6a']], site.methylations.map(&:to_a))
        assert_equal(1, blocked.fragments.size)
        assert(blocked.fragments.first.circular?)
        assert_equal(1, allowed.cuts.size)
        assert(allowed.blocked_sites.empty?)
      end

      def test_marked_pattern_defines_explicit_rule_without_database_lookup
        result = Bio::RestrictionEnzyme.digest('ggaattc', 'g^aattc')
        cut = result.cuts.first

        assert_equal(1, cut.primary_cut)
        assert_equal(5, cut.complement_cut)
        assert_equal('aatt', cut.sticky_end.sequence)
        assert_equal(2, result.fragments.size)
      end

      def test_reverse_orientation_site_is_scanned_with_reversed_methylation_map
        rule = RestrictionEnzyme::DigestRule.new('aattc',
                                                cuts: [[0, 2]],
                                                sensitive_to: [[:+, 0, nil]])
        blocked = Bio::RestrictionEnzyme.digest('tgaatt',
                                               rule,
                                               circular: true,
                                               methylation: [[:-, 0, nil]])
        cut = Bio::RestrictionEnzyme.digest('tgaatt',
                                           rule,
                                           circular: true).cuts.first

        assert(blocked.cuts.empty?)
        assert_equal(:-, blocked.blocked_sites.first.orientation)
        assert_not_nil(cut)
        assert_equal(:-, cut.orientation)
      end

      def test_reverse_orientation_maps_cuts_outside_the_recognition_site
        rule = RestrictionEnzyme::DigestRule.new('aatg', cuts: [[-1, 1]])
        result = Bio::RestrictionEnzyme.digest('acatt', rule, circular: true)

        assert_equal([[4, 2]], result.cuts.map do |cut|
          [cut.primary_cut, cut.complement_cut]
        end)
      end

      private

      def rotate(sequence, positions)
        positions += sequence.length if positions < 0
        sequence[positions..-1] + sequence[0...positions]
      end
    end
  end
end
