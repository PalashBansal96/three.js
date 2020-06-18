import {WebGLRenderTarget, WebGLRenderTargetOptions} from "./WebGLRenderTarget";
import {Texture} from "../textures/Texture";

export class WebGLMultiRenderTarget extends WebGLRenderTarget {

	constructor(
		width: number,
		height: number,
		numAttachments: number,
		options?: WebGLRenderTargetOptions
	);

	public textures: Texture[];

}
