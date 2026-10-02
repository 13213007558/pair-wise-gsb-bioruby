#
# test/unit/bio/util/restriction_enzyme/test_circular_analysis.rb -
#   Unit test for Bio::RestrictionEnzyme::CircularAnalysis
#
# Author::    BioRuby contributors
# Copyright:: Copyright (c) 2026 BioRuby project
# License::   The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 4,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/util/restriction_enzyme/circular_analysis'

# :nodoc:
module Bio
  module TestRestrictionEnzyme
    class TestCircularAnalysis < Test::Unit::TestCase # :nodoc:
      def setup
        @t = Bio::RestrictionEnzyme::CircularAnalysis
        @ecoRI = { name: 'EcoRI', site: 'g^aattc', cut_complement: 5 }
        # 'gaattc' site crosses the origin (starts at index 7)
        @plasmid = 'ttcgggggaa'
      end

      def test_linear_cut
        result = @t.cut('atgcgaattcgc', @ecoRI, circular: false)
        assert_equal({ primary: [5], complement: [9] }, result.cuts)
        assert_equal([], result.nicks)
        assert_equal([], result.blocked_sites)
        assert_equal(2, result.fragments.size)

        first, second = result.fragments
        assert_equal('atgcg    ', first.primary)
        assert_equal('tacgcttaa', first.complement)
        assert_equal(:blunt, first.left_end.type)
        assert_equal(:five_prime, first.right_end.type)
        assert_equal('aatt', first.right_end.overhang)

        assert_equal('aattcgc', second.primary)
        assert_equal('    gcg', second.complement)
        assert_equal(:five_prime, second.left_end.type)
        assert_equal('aatt', second.left_end.overhang)
        assert_equal(:blunt, second.right_end.type)
      end

      def test_circular_site_crossing_origin
        result = @t.cut(@plasmid, @ecoRI)
        assert_equal({ primary: [8], complement: [2] }, result.cuts)
        # a single double-strand cut linearizes the circle
        assert_equal(1, result.fragments.size)

        fragment = result.fragments.first
        assert_equal('aattcggggg    ', fragment.primary)
        assert_equal('    gcccccttaa', fragment.complement)
        assert_equal(:five_prime, fragment.left_end.type)
        assert_equal('aatt', fragment.left_end.overhang)
        assert_equal(:five_prime, fragment.right_end.type)
        assert_equal('aatt', fragment.right_end.overhang)
      end

      def test_circular_two_cuts_give_two_fragments
        sequence = 'gaattc' + 'cccc' + 'gaattc' + 'gg'
        result = @t.cut(sequence, @ecoRI)
        assert_equal({ primary: [1, 11], complement: [5, 15] }, result.cuts)
        assert_equal(2, result.fragments.size)
        primaries = result.fragments.map { |f| f.primary.tr(' ', '') }.sort
        assert_equal(['aattcccccg', 'aattcggg'].sort, primaries)
      end

      def test_circular_without_cuts_gives_no_fragments
        result = @t.cut('aaaaaaaa', @ecoRI)
        assert_equal([], result.fragments)
        assert_equal({ primary: [], complement: [] }, result.cuts)
      end

      def test_methylation_blocks_site
        rule = @ecoRI.merge(methylation: :blocked)
        result = @t.cut(@plasmid, rule, methylation: { primary: [8] })
        assert_equal({ primary: [], complement: [] }, result.cuts)
        assert_equal([], result.fragments)
        assert_equal(1, result.blocked_sites.size)

        site = result.blocked_sites.first
        assert_equal('EcoRI', site.enzyme)
        assert_equal(7, site.start)
        assert_equal(13, site.stop)
        assert_equal('gaattc', site.sequence)
        assert_equal([8], site.primary_methylation)
        assert_equal([], site.complement_methylation)
      end

      def test_complement_strand_methylation
        rule = @ecoRI.merge(methylation: :blocked)
        # position 0 pairs with a base inside the site (site spans 7..2)
        blocked = @t.cut(@plasmid, rule, methylation: { complement: [0] })
        assert_equal(1, blocked.blocked_sites.size)
        assert_equal([0], blocked.blocked_sites.first.complement_methylation)

        insensitive = @t.cut(@plasmid, @ecoRI, methylation: { complement: [0] })
        assert_equal([], insensitive.blocked_sites)
        assert_equal(1, insensitive.fragments.size)
      end

      def test_strand_specific_methylation_rules
        methylation = { primary: [8], complement: [0] }
        primary_only = @t.cut(@plasmid, @ecoRI.merge(methylation: :blocked_primary),
                              methylation: methylation)
        assert_equal(1, primary_only.blocked_sites.size)

        complement_only = @t.cut(@plasmid, @ecoRI.merge(methylation: :blocked_complement),
                                 methylation: methylation)
        assert_equal(1, complement_only.blocked_sites.size)

        # methylation outside the recognition site does not block
        outside = @t.cut(@plasmid, @ecoRI.merge(methylation: :blocked),
                         methylation: { primary: [4] })
        assert_equal([], outside.blocked_sites)
        assert_equal(1, outside.fragments.size)
      end

      def test_custom_methylation_callable
        # blocks only when both strands are methylated inside the site
        both = lambda do |site|
          !site.primary_methylation.empty? && !site.complement_methylation.empty?
        end
        rule = @ecoRI.merge(methylation: both)

        one_strand = @t.cut(@plasmid, rule, methylation: { primary: [8] })
        assert_equal([], one_strand.blocked_sites)
        assert_equal(1, one_strand.fragments.size)

        two_strands = @t.cut(@plasmid, rule,
                             methylation: { primary: [8], complement: [1] })
        assert_equal(1, two_strands.blocked_sites.size)
        assert_equal([], two_strands.fragments)
      end

      def test_rotation_invariance
        rule = @ecoRI.merge(methylation: :blocked)
        methylation = { primary: [8], complement: [1] }
        shift = 3
        rotated = @plasmid[shift..-1] + @plasmid[0...shift]
        rotated_methylation = {
          primary: methylation[:primary].map { |p| (p - shift) % @plasmid.size },
          complement: methylation[:complement].map { |p| (p - shift) % @plasmid.size }
        }

        # use an insensitive rule so the fragments can be compared
        insensitive = @ecoRI.merge(methylation: :none)
        original = @t.cut(@plasmid, insensitive, methylation: methylation)
        rotated_result = @t.cut(rotated, insensitive, methylation: rotated_methylation)

        n = @plasmid.size
        expected_cuts = {
          primary: original.cuts[:primary].map { |p| (p - shift) % n }.sort,
          complement: original.cuts[:complement].map { |p| (p - shift) % n }.sort
        }
        assert_equal(expected_cuts, rotated_result.cuts)

        fragment_key = lambda do |f|
          [f.primary, f.complement, f.left_end.type, f.right_end.type,
           f.left_end.overhang, f.right_end.overhang]
        end
        assert_equal(original.fragments.map(&fragment_key).sort,
                     rotated_result.fragments.map(&fragment_key).sort)

        # blocked sites are reported with rotated coordinates
        original_blocked = @t.cut(@plasmid, rule, methylation: methylation)
        rotated_blocked = @t.cut(rotated, rule, methylation: rotated_methylation)
        assert_equal(1, original_blocked.blocked_sites.size)
        assert_equal(1, rotated_blocked.blocked_sites.size)
        assert_equal((original_blocked.blocked_sites.first.start - shift) % n,
                     rotated_blocked.blocked_sites.first.start)
        assert_equal(original_blocked.blocked_sites.first.sequence,
                     rotated_blocked.blocked_sites.first.sequence)
      end

      def test_nicking_enzyme
        rule = { name: 'Nt.test', site: 'gagtc', cut_primary: 4, cut_complement: false }
        result = @t.cut('agagtca', rule, circular: false)
        assert_equal({ primary: [5], complement: [] }, result.cuts)
        assert_equal(1, result.nicks.size)
        nick = result.nicks.first
        assert_equal(:primary, nick.strand)
        assert_equal(5, nick.position)
        assert_equal('Nt.test', nick.enzyme)
        # a nick does not split a linear molecule
        assert_equal(1, result.fragments.size)
        assert_equal('agagtca', result.fragments.first.primary)
      end

      def test_complement_nicking_enzyme
        rule = { name: 'Nb.test', site: 'gagtc', cut_primary: false, cut_complement: 2 }
        result = @t.cut('agagtca', rule, circular: false)
        assert_equal({ primary: [], complement: [3] }, result.cuts)
        assert_equal(1, result.nicks.size)
        assert_equal(:complement, result.nicks.first.strand)
      end

      def test_blunt_cut_with_string_rule
        result = @t.cut('atgaattccg', 'g^aattc', circular: false)
        assert_equal({ primary: [3], complement: [3] }, result.cuts)
        assert_equal(2, result.fragments.size)
        result.fragments.each do |fragment|
          assert_equal(:blunt, fragment.left_end.type)
          assert_equal(:blunt, fragment.right_end.type)
        end
      end

      def test_three_prime_overhang
        pstI = { name: 'PstI', site: 'ctgcag', cut_primary: 5, cut_complement: 1 }
        result = @t.cut('aactgcagtt', pstI, circular: false)
        assert_equal({ primary: [7], complement: [3] }, result.cuts)
        assert_equal(2, result.fragments.size)
        assert_equal(:three_prime, result.fragments.first.right_end.type)
        assert_equal('tgca', result.fragments.first.right_end.overhang)
        assert_equal(:three_prime, result.fragments.last.left_end.type)
      end

      def test_ambiguous_site
        rule = { name: 'N-test', site: 'g^nntc', cut_complement: 3 }
        result = @t.cut('cgaatcg', rule, circular: false)
        assert_equal({ primary: [2], complement: [4] }, result.cuts)
      end

      def test_multiple_enzymes
        result = @t.cut('agaattcggagtcc', @ecoRI,
                        { name: 'blunt', site: 'gag^tc' },
                        circular: false)
        assert_equal([2, 11], result.cuts[:primary])
        assert_equal([6, 11], result.cuts[:complement])
        assert_equal(3, result.fragments.size)
      end

      def test_enzyme_rule_validation
        assert_raise(ArgumentError) { @t.normalize_rule({}) }
        assert_raise(ArgumentError) { @t.normalize_rule(site: 'gaattc') }
        assert_raise(ArgumentError) { @t.normalize_rule(42) }
      end

      def test_methylation_position_out_of_range
        assert_raise(ArgumentError) do
          @t.cut('gaattc', @ecoRI, circular: false, methylation: { primary: [6] })
        end
        # circular positions wrap around instead
        result = @t.cut('gaattc', @ecoRI, methylation: { primary: [6] })
        assert_equal(1, result.fragments.size)
      end

      def test_enzyme_rule_object
        rule = Bio::RestrictionEnzyme::CircularAnalysis::EnzymeRule.new(
          name: 'EcoRI', site: 'gaattc', cut_primary: 1, cut_complement: 5
        )
        assert_equal(1, rule.cut_primary)
        assert_equal(5, rule.cut_complement)
        assert_equal(:none, rule.methylation)
        result = @t.cut(@plasmid, rule)
        assert_equal(1, result.fragments.size)
      end
    end
  end
end
