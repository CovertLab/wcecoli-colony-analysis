from src.ridgeline import flatten


class TestFlatten:

    @staticmethod
    def test_flatten_primitive() -> None:
        flat = flatten(1)
        assert flat == [1]

    @staticmethod
    def test_flatten_str() -> None:
        flat = flatten('foo')
        assert flat == ['foo']

    @staticmethod
    def test_flatten_list() -> None:
        lst = [1, 2, 3]
        flat = flatten(lst)
        assert flat == lst

    @staticmethod
    def test_flatten_obj() -> None:
        obj = object()
        flat = flatten(obj)
        assert flat[0] is obj
        assert flat == [obj]

    @staticmethod
    def test_flatten_nested_list() -> None:
        flat = flatten([[1, 2], ['a', -1]])
        assert flat == [1, 2, 'a', -1]

    @staticmethod
    def test_flatten_depth_first() -> None:
        flat = flatten([[1, 2], [3, [4, 5]]])
        assert flat == [1, 2, 3, 4, 5]
