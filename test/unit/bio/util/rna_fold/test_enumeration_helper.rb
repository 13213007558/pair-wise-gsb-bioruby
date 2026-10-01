# frozen_string_literal: true

# test/unit/bio/util/rna_fold/test_enumeration_helper.rb - Test support load
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 4,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/util/rna_fold'
require_relative 'enumeration_reference'
