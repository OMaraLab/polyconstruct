def test_import_polyconf():
    try:
        import polyconf
    except ImportError:
        assert False, "Failed to import polyconf"
