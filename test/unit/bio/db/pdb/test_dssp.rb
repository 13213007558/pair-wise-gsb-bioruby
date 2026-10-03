# frozen_string_literal: true
#
# = test/unit/bio/db/pdb/test_dssp.rb - Unit test for Bio::PDB::DSSP
#
# Copyright:: Copyright (C) 2026 Naohisa Goto <ng@bioruby.org>
# License::   The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 4,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio'

module Bio
  # Shared helper methods for the Bio::PDB::DSSP tests.
  module PDBDSSPTestHelper
    TestDataPDB = Pathname.new(File.join(BioRubyTestDataPath,
                                         'pdb')).cleanpath.to_s

    def load_pdb(filename)
      Bio::PDB.new(File.read(File.join(TestDataPDB, filename)))
    end

    def load_dssp(filename, options = {})
      Bio::PDB::DSSP.new(load_pdb(filename), options)
    end

    # Builds a PDB ATOM record line from a spec hash with :name,
    # :resname, :chain, :resseq, :coord, :altloc and :occupancy keys.
    def atom_line(serial, spec)
      name = spec[:name].upcase
      coord = spec[:coord]
      format('ATOM  %5d %-4s%s%3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ',
             serial, name, spec[:altloc], spec[:resname], spec[:chain],
             spec[:resseq], coord[0], coord[1], coord[2], spec[:occupancy],
             20.0, name[/[A-Za-z]/])
    end

    # Builds a PDB string from residue description hashes of the form
    # { resname:, chain:, resseq:, atoms: { name => [x, y, z] or
    #   [[altloc, occupancy, x, y, z], ...] } }.
    def pdb_string(residues)
      atoms = residues.flat_map do |res|
        res[:atoms].flat_map { |name, entries| expand_entries(entries, name, res) }
      end
      atoms.each_with_index.map { |spec, i| atom_line(i + 1, spec) }.join("\n") + "\n"
    end

    def expand_entries(entries, name, res)
      entries = [[nil, nil, *entries]] if entries.first.is_a?(Numeric)
      entries.map do |altloc, occupancy, x, y, z|
        { name: name, resname: res[:resname], chain: res[:chain],
          resseq: res[:resseq], coord: [x, y, z], altloc: altloc || ' ',
          occupancy: occupancy || 1.00 }
      end
    end

    def residue(resname, chain, resseq, atoms)
      { resname: resname, chain: chain, resseq: resseq, atoms: atoms }
    end

    def dssp_from_string(str, options = {})
      Bio::PDB::DSSP.new(Bio::PDB.new(str), options)
    end
  end

  # Tests for the DSSP electrostatic hydrogen bond model and for the
  # amide hydrogen reconstruction.
  class TestPDBDSSPHydrogenBond < Test::Unit::TestCase
    include PDBDSSPTestHelper

    def test_hbond_energy_model
      oxygen = Bio::PDB::Coordinate[1.23, 0.0, 0.0]
      carbon = Bio::PDB::Coordinate[0.0, 0.0, 0.0]
      nitrogen = Bio::PDB::Coordinate[4.13, 0.0, 0.0]
      hydrogen = Bio::PDB::Coordinate[3.13, 0.0, 0.0]
      energy = Bio::PDB::DSSP.hbond_energy(oxygen, carbon, nitrogen, hydrogen)
      assert_in_delta(-2.904, energy, 0.001)
    end

    def test_hydrogen_reconstruction_geometry
      dssp = dssp_from_string(two_residue_string(3.1))
      entry = dssp.residues[1]
      n = Bio::PDB::Coordinate[3.1, 1.2, 0.2]
      prev_c = Bio::PDB::Coordinate[2.0, 1.4, 0.3]
      ca = Bio::PDB::Coordinate[4.5, 1.5, 0.3]
      bisector = ((prev_c - n).normalize + (ca - n).normalize).normalize
      expected = n - bisector * Bio::PDB::DSSP::NH_BOND_LENGTH
      assert_in_delta(1.0, (entry.h - n).r, 1.0e-6)
      3.times { |i| assert_in_delta(expected[i], entry.h[i], 1.0e-6) }
    end

    def test_first_residue_has_no_reconstructed_hydrogen
      dssp = dssp_from_string(two_residue_string(3.1))
      assert_nil(dssp.residues[0].h)
    end

    def test_no_hydrogen_reconstruction_across_chain_break
      dssp = dssp_from_string(two_residue_string(8.1))
      assert_false(dssp.residues[0].bonded_to_next)
      assert_nil(dssp.residues[1].h)
      assert_empty(dssp.residues[1].donor_hbonds)
    end

    def test_explicit_hydrogen_is_used
      atoms = { 'n' => [0.0, 0.0, 0.0], 'ca' => [1.4, 0.0, 0.0],
                'c' => [2.0, 1.3, 0.0], 'o' => [1.5, 2.4, 0.0],
                'h' => [-0.7, 0.8, 0.1] }
      dssp = dssp_from_string(pdb_string([residue('ALA', 'A', 1, atoms)]))
      h = dssp.residues[0].h
      assert_in_delta(-0.7, h[0], 1.0e-6)
      assert_in_delta(0.8, h[1], 1.0e-6)
      assert_in_delta(0.1, h[2], 1.0e-6)
    end

    def test_hbonds_between_peptide_bonded_neighbors_are_ignored
      dssp = load_dssp('model_beta_turn.pdb')
      dssp.hbonds.each do |hbond|
        assert_false(hbond.acceptor.successor.equal?(hbond.donor))
      end
    end

    def test_hbond_sorts_by_energy
      bonds = load_dssp('model_alpha_helix.pdb').hbonds.sort
      assert_operator(bonds.first.energy, :<=, bonds.last.energy)
    end

    private

    # Two alanines; the C(1)-N(2) distance depends on _n2_x_.
    def two_residue_string(n2_x)
      pdb_string([
        residue('ALA', 'A', 1, 'n' => [0.0, 0.0, 0.0], 'ca' => [1.458, 0.0, 0.0],
                'c' => [2.0, 1.4, 0.3], 'o' => [1.6, 2.5, 0.5]),
        residue('ALA', 'A', 2, 'n' => [n2_x, 1.2, 0.2], 'ca' => [n2_x + 1.4, 1.5, 0.3],
                'c' => [n2_x + 1.9, 2.9, 0.5], 'o' => [n2_x + 1.2, 3.9, 0.4])
      ])
    end
  end

  # Tests for alternate location selection, missing backbone atoms,
  # proline donors and model selection.
  class TestPDBDSSPStructureCleanup < Test::Unit::TestCase
    include PDBDSSPTestHelper

    def test_blank_altloc_is_preferred
      atoms = altloc_residue_atoms([[' ', 0.30, 1.0, 1.0, 1.0],
                                    ['A', 0.90, 9.0, 9.0, 9.0]])
      dssp = dssp_from_string(pdb_string([residue('ALA', 'A', 1, atoms)]))
      assert_in_delta(1.0, dssp.residues[0].n[0], 1.0e-6)
    end

    def test_highest_occupancy_altloc_is_used
      atoms = altloc_residue_atoms([['A', 0.40, 2.0, 2.0, 2.0],
                                    ['B', 0.60, 5.0, 5.0, 5.0]])
      dssp = dssp_from_string(pdb_string([residue('ALA', 'A', 1, atoms)]))
      assert_in_delta(5.0, dssp.residues[0].n[0], 1.0e-6)
    end

    def test_missing_oxygen_cannot_accept
      dssp = dssp_from_string(fixture_without_atom('model_beta_turn.pdb', 1, 'O'))
      assert_empty(dssp.residues[0].acceptor_hbonds)
      assert_equal('     ', dssp.secondary_structure)
    end

    def test_missing_carbonyl_carbon_breaks_hydrogen_reconstruction
      dssp = dssp_from_string(fixture_without_atom('model_beta_turn.pdb', 2, 'C'))
      assert_false(dssp.residues[1].bonded_to_next)
      assert_nil(dssp.residues[2].h)
      assert_empty(dssp.residues[2].donor_hbonds)
    end

    def test_proline_cannot_donate
      pro_str = File.read(File.join(TestDataPDB, 'model_beta_turn.pdb'))
      pro_str = pro_str.gsub('ALA A   4', 'PRO A   4')
      dssp = dssp_from_string(pro_str)
      assert_nil(dssp.residues[3].h)
      assert_empty(dssp.residues[3].donor_hbonds)
      assert_empty(dssp.residues[0].acceptor_hbonds)
      assert_equal('     ', dssp.secondary_structure)
    end

    def test_only_first_model_is_analyzed
      one = atom_line(1, atom_spec('N', 'ALA', 'A', 1, [0.0, 0.0, 0.0]))
      two = atom_line(2, atom_spec('N', 'ALA', 'A', 2, [3.0, 0.0, 0.0]))
      str = "MODEL        1\n#{one}\nENDMDL\nMODEL        2\n#{two}\nENDMDL\n"
      assert_equal(1, dssp_from_string(str).residues.size)
    end

    def test_unsupported_structure_raises
      assert_raise(ArgumentError) { Bio::PDB::DSSP.new(42) }
    end

    private

    def altloc_residue_atoms(n_entries)
      { 'n' => n_entries, 'ca' => [1.4, 0.0, 0.0],
        'c' => [2.0, 1.3, 0.0], 'o' => [1.5, 2.4, 0.0] }
    end

    def atom_spec(name, resname, chain, resseq, coord)
      { name: name, resname: resname, chain: chain, resseq: resseq,
        coord: coord, altloc: ' ', occupancy: 1.00 }
    end

    def fixture_without_atom(filename, resseq, atom_name)
      File.read(File.join(TestDataPDB, filename)).lines.reject do |line|
        line.start_with?('ATOM') && line[12..15].strip == atom_name &&
          line[22..25].to_i == resseq
      end.join
    end
  end

  # Tests for the DSSP-style assignment of the idealized model
  # fragments in test/data/pdb (alpha helix, beta strands, beta turn)
  # and for the documented conflict priority.
  class TestPDBDSSPFixtures < Test::Unit::TestCase
    include PDBDSSPTestHelper

    def setup
      @helix = load_dssp('model_alpha_helix.pdb')
      @antiparallel = load_dssp('model_beta_antiparallel.pdb')
      @parallel = load_dssp('model_beta_parallel.pdb')
      @turn = load_dssp('model_beta_turn.pdb')
    end

    def test_alpha_helix_assignment
      assert_equal(' HHHHHHHHHH ', @helix.secondary_structure)
      assert_equal(1, @helix.chains.size)
      assert_equal(12, @helix.residues.size)
    end

    def test_alpha_helix_turns_and_hbonds
      (0..7).each { |i| assert_true(@helix.residues[i].turns[4]) }
      assert_false(@helix.residues[8].turns[4])
      assert_equal(8, @helix.hbonds.size)
      @helix.hbonds.each { |hbond| assert_operator(hbond.energy, :<, -0.5) }
    end

    def test_alpha_helix_hbond_evidence
      first = @helix.residues.first
      donor = @helix.residues[4]
      assert_equal(1, first.acceptor_hbonds.size)
      assert_equal(donor, first.acceptor_hbonds.first.donor)
      assert_equal(donor.donor_hbonds.first, first.acceptor_hbonds.first)
    end

    def test_helix_takes_priority_over_turn
      entry = @helix.residues[5]
      assert_equal('H', entry.structure)
      assert_true(entry.helix?)
      assert_false(entry.turn?)
      assert_false(entry.sheet?)
    end

    def test_antiparallel_beta_assignment
      assert_equal(' EEEE  EEEE ', @antiparallel.secondary_structure)
      assert_equal(2, @antiparallel.chains.size)
    end

    def test_antiparallel_bridge_partners
      entry = @antiparallel.residues[1]
      partner, type = entry.bridges.first
      assert_equal(:antiparallel, type)
      assert_equal('B', partner.chain_id)
      assert_equal('15', partner.residue_id)
      assert_equal('E', partner.structure)
    end

    def test_ladder_takes_priority_over_bridge
      entry = @antiparallel.residues[1]
      assert_equal('E', entry.structure)
      assert_true(entry.sheet?)
      assert_false(entry.bridges.empty?)
      assert_false(entry.bridge?)
    end

    def test_parallel_beta_assignment
      assert_equal(' EEEE  EEEE ', @parallel.secondary_structure)
      partner, type = @parallel.residues[1].bridges.first
      assert_equal(:parallel, type)
      assert_equal('B', partner.chain_id)
    end

    def test_beta_turn_assignment
      assert_equal(' TT  ', @turn.secondary_structure)
      assert_true(@turn.residues[0].turns[3])
      assert_false(@turn.residues[0].turns[4])
      assert_true(@turn.residues[1].turn?)
      assert_true(@turn.residues[2].turn?)
      assert_true(@turn.residues[3].coil?)
    end

    def test_single_bridge_is_assigned_b
      dssp = dssp_from_string(single_bridge_string)
      assert_equal('B', dssp.residues[2].structure)
      assert_equal('B', dssp.residues[4].structure)
      assert_true(dssp.residues[2].bridge?)
      assert_equal(' ', dssp.residues[1].structure)
    end

    def test_single_bridge_partner
      dssp = dssp_from_string(single_bridge_string)
      partner, type = dssp.residues[2].bridges.first
      assert_equal(:antiparallel, type)
      assert_equal(dssp.residues[4], partner)
    end

    def test_chain_break_blocks_turns_and_helices
      dssp = broken_helix_dssp
      assert_false(dssp.residues[5].bonded_to_next)
      assert_nil(dssp.residues[6].h)
      assert_equal(' HHHH  HHHH ', dssp.secondary_structure)
    end

    def test_option_hbond_cutoff
      strict = load_dssp('model_alpha_helix.pdb', hbond_cutoff: -3.0)
      assert_equal(' ' * 12, strict.secondary_structure)
    end

    def test_option_max_peptide_bond
      broken = load_dssp('model_alpha_helix.pdb', max_peptide_bond: 0.1)
      assert_equal(' ' * 12, broken.secondary_structure)
      assert_empty(broken.hbonds)
    end

    def test_residue_info_accessors
      first = @helix.residues.first
      assert_equal('A', first.chain_id)
      assert_equal('ALA', first.resname)
      assert_equal('1', first.residue_id)
      assert_kind_of(Bio::PDB::Residue, first.residue)
    end

    def test_each_and_enumerable
      count = 0
      @helix.each do |info|
        assert_kind_of(Bio::PDB::DSSP::ResidueInfo, info)
        count += 1
      end
      assert_equal(12, count)
      assert_equal(12, @helix.count)
    end

    private

    # Central 3 residues of each strand of the antiparallel fixture,
    # leaving a single isolated beta bridge.
    def single_bridge_string
      chains = load_pdb('model_beta_antiparallel.pdb').models.first.chains
      keep = [[2, 3, 4], [12, 13, 14]]
      lines = chains.each_with_index.flat_map do |chain, ci|
        chain.residues.select { |res| keep[ci].include?(res.resSeq) }
             .flat_map { |res| res.atoms.map(&:to_s) }
      end
      lines.join
    end

    # Alpha helix fixture with the nitrogen of residue 7 moved 20 A
    # away, breaking the peptide chain between residues 6 and 7.
    def broken_helix_dssp
      residues = load_pdb('model_alpha_helix.pdb').models.first.chains.first.residues
      nitrogen = residues[6].atoms.find { |atom| atom.name == 'N' }
      nitrogen.x = nitrogen.x + 20.0
      Bio::PDB::DSSP.new(residues)
    end
  end
end
