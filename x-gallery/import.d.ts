import { InlineHandler } from './handler';
export declare const promise_all: () => Promise<any[]>;
export declare class ImportHandler implements InlineHandler {
    pattern: RegExp;
    make_span: (body_raw: string) => any;
}
