import { InlineHandler } from './handler';
export declare class InlineRefHandler implements InlineHandler {
    pattern: RegExp;
    make_span: (body_raw: string) => HTMLSpanElement;
}
