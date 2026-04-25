from __future__ import annotations

import os

DEFAULT_HOST = "Escherichia coli"
DEFAULT_TAG_TYPE = "6xHis"
DEFAULT_PROTEASE_SITE = "Enterokinase"
DEFAULT_VECTOR = "pET-28a(+)"
DEFAULT_CLONING_SITES = ("BamHI", "XhoI")
DEFAULT_MAX_REMEDIATION_ATTEMPTS = 3
LLM_MODEL = os.getenv("LLM_MODEL", "deepseek/deepseek-v4-flash")
LLM_MAX_TOKENS = int(os.getenv("LLM_MAX_TOKENS", "4096"))
OPENROUTER_BASE_URL = "https://openrouter.ai/api/v1"
OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY", "")

from langchain_openai import ChatOpenAI

class ChatOpenRouter(ChatOpenAI):
    """
    Wrapper around ChatOpenAI that ensures max_tokens is passed instead of max_completion_tokens.
    OpenRouter uses max_tokens to pre-calculate credit limits. If it encounters max_completion_tokens
    (which LangChain automatically translates to for certain models), OpenRouter assumes the maximum
    context length (e.g. 65536) and may reject the request with a 402 Error.
    """
    def _get_request_payload(self, messages, *, stop = None, **kwargs):
        payload = super()._get_request_payload(messages, stop=stop, **kwargs)
        if "max_completion_tokens" in payload:
            payload["max_tokens"] = payload.pop("max_completion_tokens")
        return payload
