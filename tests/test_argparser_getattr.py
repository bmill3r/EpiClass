import argparse


class ParserClass(object):
    def __init__(self, input_args=None):

        """INITIALIZE PARSER"""
        parser = argparse.ArgumentParser()
        parser.add_argument('--name', type=str, required=False, default='World')
        parser.add_argument(
            '--testFlag',
            action='store_true',
            required=False)
        self.args = parser.parse_args(input_args)
        print(self.args)

    def run(self):
        # self.args is the arg Namespace from parser.
        # for each arg, if it's boolean and True, then use getattr to run the method whose name matches the arg key (name of arg).
        return [getattr(self, key)() for key, val in vars(self.args).items() if isinstance(val, bool) and val]


    def testFlag(self):
        return "Hello, {}!".format(self.args.name)


def test_parser_getattr():
    assert ["Hello, EpiClass!"] == ParserClass("--testFlag --name EpiClass".split()).run()


if __name__ == '__main__':
    ParserClass().run()
